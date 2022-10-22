#!/usr/bin/env perl
use strict;
# use warnings;
use 5.010;

# 用户输入待编辑的染色体($chr)和起始位点($start_position)以及编辑模式($edit_pattern),其中编辑模式以/或+或-三种符号开头，后面是具体的碱基序列。
# 若编辑模式以/开头，代表碱基替换，例如/A，则表示将染色体起始位点上的碱基替换成A（被替换的碱基和替换的碱基不能相同，否则操作无意义，后面的代码会进行检查和提示）
# 若编辑模式以+开头，代表碱基插入，例如+ATGCC，则表示在染色体起始位点的碱基后插入ATGCC
# 若编辑模式以-开头，代表碱基删除，例如-GAG，则表示在染色体起始位点的碱基后删除GAG（被删除的序列需要本来就存在，不然无法进行删除，后面的代码会进行检查和提示）
# 例如：对于用户输入的染色体起始位点chr1-643995（643995是碱基t，其后的连续20个碱基分别是ctccttctagagacatggta）
# 如果要将chr1-643995上的t替换成A,那么用户输入的编辑模式为/A；
# 如果要在chr1-643995上的t之后插入GTATT,那么用户输入的编辑模式为+GTATT；
# 如果要将chr1-643995上的t之后的ctccttc删除,那么用户输入的编辑模式为-ctccttc；

# 另外对于PAM($pam_type)，让用户可以在网页选择是NGG还是NG这两种的哪种(默认为NGG)
# 对于切口位点到编辑位点的最大距离($dis_nickase)，让用户可以在网页选择设定（默认为10）
# 对于PBS的长度范围（$min_len_PBS，$max_len_PBS），让用户可以在网页选择设定（默认分别为8，17）
# 对于RT的长度范围（$min_len_RT，$max_len_RT），让用户可以在网页选择设定（默认分别为8，24）
# 对于sgRNA的距离范围（$min_dis_sgRNA_nickase,$max_dis_sgRNA_nickase），让用户可以在网页选择设定（默认分别为0，100）
# 对于RT下游的最短同源序列长度$homology，让用户可以在网页选择设定（默认5）

##以下几个参数为用户输入的，否则采用下面的默认参数
#my $chr='chr1';	#染色体
#my $start_position=643995;	#染色体上的起始位置
#my $edit_pattern='-ctccttc';	#编辑模式, /T表示将当前碱基替换成T, +GTATT表示在当前碱基后插入GTATT，-ctccttc表示在当前碱基后删除子序列ctccttc
#
#my $pam_type='NGG';	#PAM	选择实验用的是3bp的NGG作为PAM，N是任意还是2bp的NG作为PAM。
#my $dis_nickase=10;	#maximum edit-to-nick distance 切口位点到编辑位点的最大距离，根据实验设置在10bp以内，但是发现有一半的位点在10bp内没有PAM，故改成50bp
#my ($min_dis_sgRNA_nickase,$max_dis_sgRNA_nickase)=(0,100);	#minimum/maximum nick-to-nick distance	PE3的第二个sgRNA切口位点到pegRNA切口位点的最大距离，根据David Liu和PrimeDesign论文设置为0~100
#my ($max_len_PBS,$max_len_RT)=(17,24);	#实验设置的RT模板和PBS最长长度
#my ($min_len_PBS,$min_len_RT)=(8,8);	#实验设置的RT模板和PBS最短长度
#my $homology=5;	#minimum downstream homology 	RT最右端至少要与编辑位点最右侧有5bp的距离（如果距离太近是没有功能的）
my $chr=$ARGV[0];
my $start_position=$ARGV[1];
my $edit_pattern=$ARGV[2];
my $pam_type=$ARGV[3];
my $dis_nickase=$ARGV[4];
my $min_dis_sgRNA_nickase=$ARGV[5];
my $max_dis_sgRNA_nickase=$ARGV[6];
my $max_len_PBS=$ARGV[7];
my $max_len_RT=$ARGV[8];
my $min_len_PBS=$ARGV[9];
my $min_len_RT=$ARGV[10];
my $homology=$ARGV[11];


##此部分为固定的参数
my $type='User';
my $pam=substr($pam_type,1);	#去掉PAM中的N，便于后面匹配
my $length=100;
my $len_spacer=20;	#spacer长度固定20bp
my $len_target=47;	#Total 47 bps = 4 bp neighboring sequence + 20 bp protospacer + 3 bp NGG PAM+ 20 bp neighboring sequence
my $len_up=4;	#NBT的spacer上游序列的长度
my $dir_out="Position";
mkdir $dir_out unless -e $dir_out;
open(OUT,">","Position/pegRNA_April_pos.$type.txt") or die($!);
say OUT "PegRNAID\tStrand\tTarget\tSpacer\tPAM\tPBS\tRT\tEditToNickDistance\tsgRNASpacer\tsgRNAPAM\tNickToNickDistance";

my $seq=extract_sequence($chr,$start_position,$length);
die("The reference sequence of $chr:$start_position does not exist!!!\n") unless defined($seq);
my $user_seq=$seq;	#用户序列
if($edit_pattern=~/^\/[acgtACTG]$/)
{
	substr($user_seq,$length+1,0)="($edit_pattern)";
}
elsif($edit_pattern=~/^\+[acgtACTG]+$/)
{
	substr($user_seq,$length+1,0)="($edit_pattern)";
}
elsif($edit_pattern=~/^-[acgtACTG]+$/)
{
	my $temp=$edit_pattern;
	$temp=~s/^-//;
	die("The following reference sequence after $chr:$start_position is ".substr($user_seq,$length+1,20).", but the deleted sequence $temp does not exist and thus can not be deleted!!\n") unless uc(substr($user_seq,$length+1,length($temp))) eq uc($temp);
	substr($user_seq,$length+1,length($temp))="($edit_pattern)";
}
# say "The final sequence from $chr $start_position for $edit_pattern is: \n$user_seq";

my $user_seq_minus='';		#负链
my ($up, $mid, $down)=("","","");
die("The flanking sequence should only contain character of 'acgtACTG', but contain character such as N!!!\n") unless $user_seq=~/^[acgtACTG+-\/()]+$/;
if($user_seq=~/^[acgtACTG]*\(-[acgtACTG]+\)[acgtACTG]*$/)	#对于包含-的删除操作
{
	($up, $mid, $down) = $user_seq=~/^([acgtACTG]*)\(-([acgtACTG]+)\)([acgtACTG]*)$/;
	$user_seq_minus=rev(complement($down)).'(-'.rev(complement($mid)).')'.rev(complement($up));
}
elsif($user_seq=~/^[acgtACTG]*\(\+[acgtACTG]+\)[acgtACTG]*$/)	#对于包含+的插入操作
{
	($up, $mid, $down) = $user_seq=~/^([acgtACTG]*)\(\+([acgtACTG]+)\)([acgtACTG]*)$/;
	$user_seq_minus=rev(complement($down)).'(+'.rev(complement($mid)).')'.rev(complement($up));
}
elsif($user_seq=~/^[acgtACTG]+\(\/[acgtACTG]\)[acgtACTG]*$/)	#对于包含/的替换操作
{
	($up, $mid, $down) = $user_seq=~/^([acgtACTG]+)\(\/([acgtACTG])\)([acgtACTG]*)$/;
	my $seq_up=rev(complement($up));
	my $s=substr($seq_up,0,1);
	$seq_up=substr($seq_up,1);
	$user_seq_minus=rev(complement($down)).$s.'(/'.rev(complement($mid)).')'.$seq_up;
}
else
{
	my $example='
	Example:
	Input position: chr1:643995
	Sequence started at the input position is tctccttctagagacatggta
	Example for substitution: /A to substitute the base of the input position with A
	Example for inserttion: +ATT to insert ATT after the base of the input position
	Example for deletiion: -ctccttc to delete ctccttc after the base of the input position
	';
	die("The format of the intended edit is incorrect!! Please check and input a correct one!!\n$example\n");
}
# die('There should be at least one nucleotide before the parentheses') if($up eq "");
# die('There should be at least one nucleotide after the parentheses') if($down eq "");

my $count_variants=0;
my $count_rnas=0;

my @pegrnas=search_candicante($user_seq,0,'+');	#正链
if(@pegrnas)
{
	$count_rnas+=@pegrnas;
	foreach my $rna(@pegrnas)
	{
		say OUT "$rna";
	}
}

my $id=0;
$id=(split(/\t/,$pegrnas[-1]))[0] if @pegrnas;	#获取正链上的最后一个pegRNA的ID
my @pegrnas_minus=search_candicante($user_seq_minus,$id,'-');	#负链
if(@pegrnas_minus)
{
	$count_rnas+=@pegrnas_minus;
	foreach my $rna(@pegrnas_minus)
	{
		say OUT "$rna";
	}
}
$count_variants++ if @pegrnas or @pegrnas_minus;
say "$count_rnas pegRNAs";



sub search_candicante
{
	my $user_seq=$_[0];
	my $id_peg=$_[1];	#peg的起始id，默认为0
	my $strand=$_[2];	#正负链的标记，+-
	my $ref_seq=$user_seq;
	if($user_seq=~/^[acgtACTG+-\/()]+$/)	#检查，输入序列只能包含actgATCG以及()+-/
	{
		my $edit_pos =index($user_seq,'(')-1;	# 0 for the first bp
		# die("The current length of upstream sequence is ".($edit_pos+1)."bp and should be at least ".($len_spacer-3)."bp.") if $edit_pos+1<=$len_spacer-3;	#编辑位点上游至少要有17bp的长度，17是下界，最好是更长。
		
		my ($reference_allele, $alternate_allele, $type)=("","","");
		if($user_seq=~/\(-[acgtACTG]+\)/)	#对于包含-的删除操作,检查是否按合法的格式
		{
			($reference_allele) = $user_seq=~/[acgtACTG]\(-([acgtACTG]+)\)/;
			$type='-';
			$ref_seq=~s/[^acgtACTG]//g;
		}
		elsif($user_seq=~/\(\+[acgtACTG]+\)/)	#对于包含+的插入操作,检查是否按合法的格式
		{
			($alternate_allele) = $user_seq=~/[acgtACTG]\(\+([acgtACTG]+)\)/;
			$type='+';
			$ref_seq=~s/\(.+\)//g;
		}
		elsif($user_seq=~/\(\/[acgtACTG]\)/)	#对于包含/的替换操作,检查是否按合法的格式
		{
			($reference_allele, $alternate_allele) = $user_seq=~/([acgtACTG])\(\/([acgtACTG])\)/;
			$type='/';
			if(uc($reference_allele) ne uc($alternate_allele))	#对于包含/的替换操作,替换的碱基需与替换前的碱基不同，不然就不需要做任何
			{
				$ref_seq=~s/\(.+\)//g;
			}
			else
			{
				die("The substitution nucleotide is the same!!!\n");
			}
		}
		else
		{
			die("The format of input sequence is not correct!!!\n");
		}
		
		$edit_pos++ unless $type eq '/';
		my $alt_seq= $ref_seq;
		substr($alt_seq,$edit_pos,length($reference_allele))=$alternate_allele;
		
		# say "Input sequence:";
		# say $user_seq;
		# say "Original sequence:";
		# say $ref_seq;
		# say "Edited sequence:";
		# say $alt_seq;
		
		my $change=0;
		$change=abs(length($alternate_allele)-length($reference_allele)) if $type eq '+';	
		my $sub_seq=uc(substr($ref_seq,0,$edit_pos+4+length($pam)));
		my $i=rindex($sub_seq,$pam);	#NGG PAM中点的坐标或者 NG PAM中G的坐标
		my @locs;	#pegRNA的NGG PAM中点的坐标或者 NG PAM中G的坐标
		while( $i>=$len_spacer+1 and $edit_pos-($i-5)<=$dis_nickase )	#$i>=$len_spacer+1保证在PAM前有20bp的spacer；$i-5是切口位点位置的下标，$edit_pos是编辑位点的下标
		{
			push @locs,$i;
			$i=rindex($sub_seq,$pam,$i-1);		
		}
		
		my $sub_seq_pe3=uc(substr($ref_seq,0,length($ref_seq)-($len_spacer+1)));
		my $i_pe3=index($sub_seq_pe3,rev(complement($pam)));	#CCN或CN中第一个C的坐标
		my @locs_pe3;	#PE3 sgRNA的CCN或CN中第一个C的坐标
		while( $i_pe3>=0 )
		{
			push @locs_pe3,$i_pe3;
			$i_pe3=index($sub_seq_pe3,rev(complement($pam)),$i_pe3+1);
		}
		
		my @pegrnas;
		# die("There is no PAM and spacer within ${dis_nickase}bp from the edit site!!!") unless @locs;	#没有找到对应PAM
		return @pegrnas unless @locs;	#没有找到对应PAM就跳过
		foreach my $i_pam(@locs)
		{
			my $best_sg=-1;	#
			my @locs_pe3_filtered;	#满足距离条件的PE3 sgRNA的CCN或CN中第一个C的坐标
			foreach my $i_pam_sg(@locs_pe3)
			{
				my $nick_to_peg_dis=($i_pam_sg+3+length($pam))-($i_pam-5);
				next if $nick_to_peg_dis<$min_dis_sgRNA_nickase or $nick_to_peg_dis>$max_dis_sgRNA_nickase;
				
				push @locs_pe3_filtered,$i_pam_sg;
				my $ind_sg_l=$i_pam_sg+length($pam)+1;	#sgRNA 20nt spacer最左侧的碱基的绝对坐标
				my $ind_sg_r=$i_pam_sg+length($pam)+$len_spacer;	#sgRNA 20nt spacer最右侧的碱基的绝对坐标
				if($ind_sg_l<=$edit_pos and $edit_pos<=$ind_sg_r)	#PE3b
				{
					$best_sg=$i_pam_sg;
					# last;
				}
				
			}
			my $seq_sg;
			my $seq_sg_pam;
			my $nick_to_peg_dis;
			if($best_sg>=0)	#PE3b
			{
				$seq_sg=substr($alt_seq,$best_sg,$len_spacer+length($pam)+1);
				$seq_sg=rev(complement($seq_sg));
				$seq_sg_pam=substr($seq_sg,-length($pam)-1);
				$seq_sg=substr($seq_sg,0,$len_spacer);
				$nick_to_peg_dis=($best_sg+3+length($pam))-($i_pam-5);
			}	
			elsif($best_sg<0 and @locs_pe3_filtered>0)	#PE3
			{
				@locs_pe3_filtered=sort { abs(($a+3+length($pam))-($i_pam-5)-75) <=> abs(($b+3+length($pam))-($i_pam-5)-75) } @locs_pe3_filtered;
				$best_sg=$locs_pe3_filtered[0];	#根据primedesign这篇文章的排序规则，PE3 annotation at a distance as close to 75 bp away
				$seq_sg=substr($ref_seq,$best_sg,$len_spacer+length($pam)+1);
				$seq_sg=rev(complement($seq_sg));
				$seq_sg_pam=substr($seq_sg,-length($pam)-1);
				$seq_sg=substr($seq_sg,0,$len_spacer);
				$nick_to_peg_dis=($best_sg+3+length($pam))-($i_pam-5);
			}
			else	# no sgRNA
			{
				$seq_sg='Na';
				$seq_sg_pam='Na';
				$nick_to_peg_dis='Na';
			}
			
			my $pos=$edit_pos-($i_pam-5)+$change;	#pos是替换序列最右侧碱基的坐标,$edit_pos-($i_pam-5)是编辑位点相对切口位点的坐标，一般都是+1，+2，...
			$pos-- unless $type eq '/';
			foreach my $len_PBS($min_len_PBS..$max_len_PBS)
			{
				foreach my $len_RT($min_len_RT..$max_len_RT)
				{
					next unless $len_RT>=$pos+$homology;	#RT最右端至少要与编辑位点最右侧有5bp的距离（如果距离太近是没有功能的）

					my $seq_spacer=substr($ref_seq,$i_pam-1-$len_spacer,$len_spacer);
					my $seq_pam=substr($ref_seq,$i_pam-1,length($pam)+1);
					my $seq_PBS=substr($ref_seq,$i_pam-4-$len_PBS,$len_PBS);
					my $seq_RT=substr($alt_seq,$i_pam-4,$len_RT);
					my $seq_target;
					if($i_pam-1-$len_spacer-$len_up>=0)
					{
						$seq_target=substr($ref_seq,$i_pam-1-$len_spacer-$len_up,$len_target);	#这里的seq_target序列为Spacer上游序列+Spacer+PAM+Spacer下游游序列，其中Total 47 bps = 4 bp neighboring sequence + 20 bp protospacer + 3 bp NGG PAM+ 20 bp neighboring sequence
					}
					else
					{
						$seq_target=substr($ref_seq,0,$len_target);
					}
					
					$seq_PBS=complement($seq_PBS);	#PBS与DNA正链互补
					$seq_PBS=rev($seq_PBS);	#PBS调整为5'端开始
					$seq_RT=complement($seq_RT);	#RT与DNA正链互补
					$seq_RT=rev($seq_RT);	#RT调整为5'端开始
					
					$id_peg++;
					push @pegrnas,"$id_peg\t$strand\t$seq_target\t$seq_spacer\t$seq_pam\t$seq_PBS\t$seq_RT\t".($edit_pos-($i_pam-5))."\t$seq_sg\t$seq_sg_pam\t$nick_to_peg_dis";
					# say OUT "$id_peg\t$seq_target\t$seq_spacer\t$seq_PBS\t$seq_RT\t+".($edit_pos-($i_pam-5));
				}
			}
		}
		return @pegrnas;
		
	}
	else
	{
		die("The input sequence should only contain character of 'acgtACTG+-/()'!!!\n");
	}
}


sub complement	
{
	#取互补链
	my $s=$_[0];
	my %maps=('A','T','T','A','C','G','G','C');
	my @temp;
	foreach (split(//,$s))
	{
		push @temp,$maps{uc($_)};
	}
	join('',@temp);
}


sub rev	
{
	#取反向链
	join('',reverse(split(//,$_[0])));
}


sub extract_sequence
{
	my $chr;
	my $mid;
	my $len;
	if(@_==3)
	{
		$chr=$_[0];
		$mid=$_[1];
		$len=$_[2];
	}
	else
	{
		$chr='chr1';
		$mid=1014143;
		$len=2;
	}
	my ($start,$end)=($mid-$len,$mid+$len);
	# open(FILE_OUT,">","chr_location.fa");

	# print "$chr:$start:$end\n";

	my $num_colunm=50;
	open(FILE,"<","././Reference_Genome/$chr.fa") or die("The input chromosome is not included, and please change to another chromosome!\n");
	my $row1=int(($start-1)/$num_colunm)+1;
	my $row2=int(($end-1)/$num_colunm)+1;
	my $length=$end-$start+1;
	my $i=$row1;
	my $data;
	seek(FILE,0,0);
	<FILE>;#skip the first line
	seek(FILE,$start+($row1-2),1); #-1 for 0-based bed  or  -2 for 1-based bed
	while(<FILE>){
		  s/[\r\n]+$//;
		  if($row1==$row2){
			  $data.=substr($_,0,$length);
			  last;
		  }
		  elsif($i<$row2){
			   $data.=$_;
			   $i++;
		  }
		  else{
			   $data.=substr($_,0,($end-1)%$num_colunm+1);
			   last;
		  }
	}
	# print FILE_OUT ">$chr-$start-$end\n$data\n";
	# print ">$chr\t$start\t$end\n$data\n";
	close FILE;
	return $data;
}