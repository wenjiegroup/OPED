#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;


# Example for substitution: CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTGGCGCGA(/T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC
# Example for inserttion: CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTG(+ATT)GCGCGAGGCCGCCTGGCAACTCTGCGACTACTACCTGCC
# Example for deletiion: CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGC(-GCTGGCGCGA)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC
# 如上的三个例子代表了三种不同编辑类型的输入格式，三种编辑类型需要编辑的碱基都是用一对圆括号()标注。
# 对于碱基替换则用反斜杠/表示,例子表示将A替换成T；
# 对于碱基插入则用加号+表示,例子表示将在圆括号位置插入ATT；
# 对于碱基删除则用加号-表示,例子表示将在圆括号位置删除GCTGGCGCGA；

# 用户输入的序列($user_seq)需满足上面说的格式
# 另外对于PAM($pam_type)，让用户可以在网页选择是NGG还是NG这两种的哪种(默认为NGG)
# 对于切口位点到编辑位点的最大距离($dis_nickase)，让用户可以在网页选择设定（默认为10）
# 对于PBS的长度范围（$min_len_PBS，$max_len_PBS），让用户可以在网页选择设定（默认分别为8，17）
# 对于RT的长度范围（$min_len_RT，$max_len_RT），让用户可以在网页选择设定（默认分别为8，24）



my $user_seq=ARGV[0];
my $pam_type=ARGV[1];
my $dis_nickase=ARGV[2];
my $max_len_PBS=ARGV[3];
my $max_len_RT=ARGV[4];
my $min_len_PBS=ARGV[5];
my $min_len_RT=ARGV[6];

##此部分为固定的参数
my $type='User';
my $pam=substr($pam_type,1);	#去掉PAM中的N，便于后面匹配
my $len_spacer=20;	#spacer长度固定20bp
my $len_target=47;	#Total 47 bps = 4 bp neighboring sequence + 20 bp protospacer + 3 bp NGG PAM+ 20 bp neighboring sequence
my $len_up=4;	#NBT的spacer上游序列的长度


###以上师兄新版本

##以下几个参数为用户输入的，否则采用下面的默认参数
#my $user_seq='CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTGGCGCGA(/T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC';	#用户输入序列
#
#my $pam_type='NGG';	#选择实验用的是3bp的NGG作为PAM，N是任意还是2bp的NG作为PAM。
#my $dis_nickase=10;	#切口位点到编辑位点的最大距离，根据实验设置在10bp以内，但是发现有一半的位点在10bp内没有PAM，故改成50bp
#my ($max_len_PBS,$max_len_RT)=(17,24);	#实验设置的RT模板和PBS最长长度
#my ($min_len_PBS,$min_len_RT)=(8,8);	#实验设置的RT模板和PBS最短长度

##此部分为固定的参数
my $type='User';
my $pam=substr($pam_type,1);	#去掉PAM中的N，便于后面匹配
my $len_spacer=20;	#spacer长度固定20bp
my $len_target=47;	#Total 47 bps = 4 bp neighboring sequence + 20 bp protospacer + 3 bp NGG PAM+ 20 bp neighboring sequence
my $len_up=4;	#NBT的spacer上游序列的长度
my $dir_out="Sequence";
mkdir $dir_out unless -e $dir_out;
open(OUT,">","$dir_out/pegRNA.$type.txt") or die($!);
say OUT "PegRNAID\tTarget(47bp)\tSpacer\tPBS\tRT\tNickingCoordinate";



my $ref_seq=$user_seq;
if($user_seq=~/^[acgtACTG+-\/()]+$/)	#检查，输入序列只能包含actgATCG以及()+-/
{
	my $edit_pos =index($user_seq,'(')-1;	# 0 for the first bp
	die("The current length of upstream sequence is ".($edit_pos+1)."bp and should be at least ".($len_spacer-3)."bp.") if $edit_pos+1<=$len_spacer-3;	#编辑位点上游至少要有17bp的长度，17是下界，最好是更长。
	
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

		     @HJS
		    say OUT "wrongsame\twrongsame\twrongsame\twrongsame\twrongsame\twrongsome"
			die("The substitution nucleotide is the same!!!");
#

		}
	}
	else
	{
		die("The format of input sequence is not correct!!!");
	}
	
	$edit_pos++ unless $type eq '/';
	my $alt_seq= $ref_seq;
	substr($alt_seq,$edit_pos,length($reference_allele))=$alternate_allele;
	
	say "Input sequence:";
	say $user_seq;
	say "Original sequence:";
	say $ref_seq;
	say "Edited sequence:";
	say $alt_seq;
	
	my $change=0;
	$change=abs(length($alternate_allele)-length($reference_allele)) if $type eq '+';	
	my $sub_seq=uc(substr($ref_seq,0,$edit_pos+4+length($pam)));
	my $i=rindex($sub_seq,$pam);	#NGG PAM中点的坐标
	my @locs;	
	while( $i>=$len_spacer+1 and $edit_pos-($i-5)<=$dis_nickase )	#$i>=$len_spacer+1保证在PAM前有20bp的spacer；$i-5是切口位点位置的下标，$edit_pos是编辑位点的下标
	{
		push @locs,$i;
		$i=rindex($sub_seq,$pam,$i-1);		
	}
	die("There is no PAM and spacer within ${dis_nickase}bp from the edit site!!!") unless @locs;	#没有找到对应PAM
	my $id_peg=0;
	foreach my $i_pam(@locs)
	{
		my $pos=$edit_pos-($i_pam-5)+$change;	#pos是替换序列最右侧碱基的坐标,$edit_pos-($i_pam-5)是编辑位点相对切口位点的坐标，一般都是+1，+2，...
		foreach my $len_PBS($min_len_PBS..$max_len_PBS)
		{
			foreach my $len_RT($min_len_RT..$max_len_RT)
			{
				next unless $len_RT>=$pos+5;	#RT最右端至少要与编辑位点最右侧有5bp的距离（如果距离太近是没有功能的）

				my $seq_spacer=substr($ref_seq,$i_pam-1-$len_spacer,$len_spacer);
				my $seq_PBS=substr($ref_seq,$i_pam-4-$len_PBS,$len_PBS);
				my $seq_RT=substr($alt_seq,$i_pam-4,$len_RT);
				my $seq_target=substr($ref_seq,$i_pam-1-$len_spacer-$len_up,$len_target);	#这里的seq_target序列为Spacer上游序列+Spacer+PAM+Spacer下游游序列，其中Total 47 bps = 4 bp neighboring sequence + 20 bp protospacer + 3 bp NGG PAM+ 20 bp neighboring sequence
				
				$seq_PBS=complement($seq_PBS);	#PBS与DNA正链互补
				$seq_PBS=rev($seq_PBS);	#PBS调整为5'端开始
				$seq_RT=complement($seq_RT);	#RT与DNA正链互补
				$seq_RT=rev($seq_RT);	#RT调整为5'端开始
				
				$id_peg++;
				say OUT "$id_peg\t$seq_target\t$seq_spacer\t$seq_PBS\t$seq_RT\t+".($edit_pos-($i_pam-5));
			}
		}
	}
	
	
}
else
{
	die("$!: The input sequence should only contain character of 'acgtACTG+-/()'!!!");
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