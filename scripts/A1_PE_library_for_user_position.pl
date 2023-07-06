#!/usr/bin/env perl
use strict;
# use warnings;
use 5.010;

#  xuhao =====================
# my $Chromosome='chr1';
# my $Position=643995;
# my $Pattern='-ctccttc';
# my $PAM='NGG';
# my $CUT_SIZE=10;	#maximum edit-to-nick distance 
# my ($MIN_DISGRNA,$MAX_DISGRNA)=(0,100);	#minimum/maximum nick-to-nick distance
# my ($MAX_PBS,$MAX_RT)=(17,24);
# my ($MIN_PBS,$MIN_RT)=(8,8);
# my $HOMOLOGY=5;	#minimum downstream HOMOLOGY
# chdir '../';
my $Chromosome=$ARGV[0];
my $Position=$ARGV[1];
my $Pattern=$ARGV[2];
my $PAM=$ARGV[3];
my $CUT_SIZE=$ARGV[4];
my $MIN_DISGRNA=$ARGV[5];
my $MAX_DISGRNA=$ARGV[6];
my $MAX_PBS=$ARGV[7];
my $MAX_RT=$ARGV[8];
my $MIN_PBS=$ARGV[9];
my $MIN_RT=$ARGV[10];
my $HOMOLOGY=$ARGV[11];
my $time=$ARGV[12];
#  xuhao =====================


my $type='User';
my $pam=substr($PAM,1);
my $length=100;
my $len_spacer=20;	#spacer 20bp
my $len_target=47;	#Total 47 bps = 4 bp neighboring sequence + 20 bp protospacer + 3 bp NGG PAM+ 20 bp neighboring sequence
my $len_up=4;	#NBT spacer up stream
my $dir_out="Temp";
mkdir $dir_out unless -e $dir_out;
open(OUT,">","$dir_out/Position.request.$type.$time.txt") or die($!);
say OUT "PegRNAID\tStrand\tTarget\tSpacer\tPAM\tPBS\tRT\tEditToNickDistance\tsgRNASpacer\tsgRNAPAM\tNickToNickDistance";

my $seq=extract_sequence($Chromosome,$Position,$length);
die("The reference sequence of $Chromosome:$Position does not exist!!!\n") unless defined($seq);
my $SEQUENCE=$seq;
if($Pattern=~/^\/[acgtACTG]$/)
{
	substr($SEQUENCE,$length+1,0)="($Pattern)";
}
elsif($Pattern=~/^\+[acgtACTG]+$/)
{
	substr($SEQUENCE,$length+1,0)="($Pattern)";
}
elsif($Pattern=~/^-[acgtACTG]+$/)
{
	my $temp=$Pattern;
	$temp=~s/^-//;
	die("The following reference sequence after $Chromosome:$Position is ".substr($SEQUENCE,$length+1,20).", but the deleted sequence $temp does not exist and thus can not be deleted!!\n") unless uc(substr($SEQUENCE,$length+1,length($temp))) eq uc($temp);
	substr($SEQUENCE,$length+1,length($temp))="($Pattern)";
}
# say "The final sequence from $Chromosome $Position for $Pattern is: \n$SEQUENCE";

my $user_seq_minus='';		# - Strand
my ($up, $mid, $down)=("","","");
die("The flanking sequence should only contain character of 'acgtACTG', but contain character such as N!!!\n") unless $SEQUENCE=~/^[acgtACTG+-\/()]+$/;
if($SEQUENCE=~/^[acgtACTG]*\(-[acgtACTG]+\)[acgtACTG]*$/)	# - 
{
	($up, $mid, $down) = $SEQUENCE=~/^([acgtACTG]*)\(-([acgtACTG]+)\)([acgtACTG]*)$/;
	$user_seq_minus=rev(complement($down)).'(-'.rev(complement($mid)).')'.rev(complement($up));
}
elsif($SEQUENCE=~/^[acgtACTG]*\(\+[acgtACTG]+\)[acgtACTG]*$/)	# + 
{
	($up, $mid, $down) = $SEQUENCE=~/^([acgtACTG]*)\(\+([acgtACTG]+)\)([acgtACTG]*)$/;
	$user_seq_minus=rev(complement($down)).'(+'.rev(complement($mid)).')'.rev(complement($up));
}
elsif($SEQUENCE=~/^[acgtACTG]+\(\/[acgtACTG]\)[acgtACTG]*$/)	# / 
{
	($up, $mid, $down) = $SEQUENCE=~/^([acgtACTG]+)\(\/([acgtACTG])\)([acgtACTG]*)$/;
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
	die("The format of the intended edit is incorrect!! Please check and input a correct one!!\n");
}
# die('There should be at least one nucleotide before the parentheses') if($up eq "");
# die('There should be at least one nucleotide after the parentheses') if($down eq "");

my $count_variants=0;
my $count_rnas=0;

my @pegrnas=search_candicante($SEQUENCE,0,'+');	# + strand
if(@pegrnas)
{
	$count_rnas+=@pegrnas;
	foreach my $rna(@pegrnas)
	{
		say OUT "$rna";
	}
}

my $id=0;
$id=(split(/\t/,$pegrnas[-1]))[0] if @pegrnas;	# pegRNA ID
my @pegrnas_minus=search_candicante($user_seq_minus,$id,'-');	# - stand
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
	my $SEQUENCE=$_[0];
	my $id_peg=$_[1];	#peg start id start from 0
	my $Strand=$_[2];	#stand: +-
	my $ref_seq=$SEQUENCE;
	if($SEQUENCE=~/^[acgtACTG+-\/()]+$/)
	{
		my $edit_pos =index($SEQUENCE,'(')-1;	# 0 for the first bp
		die("The current length of upstream sequence is ".($edit_pos+1)."bp and should be at least ".($len_spacer-3)."bp.") if $edit_pos+1<=$len_spacer-3;	# at least 17bp upper
		
		my ($reference_allele, $alternate_allele, $type)=("","","");
		if($SEQUENCE=~/\(-[acgtACTG]+\)/)	# - 
		{
			($reference_allele) = $SEQUENCE=~/[acgtACTG]\(-([acgtACTG]+)\)/;
			$type='-';
			$ref_seq=~s/[^acgtACTG]//g;
		}
		elsif($SEQUENCE=~/\(\+[acgtACTG]+\)/)	# + 
		{
			($alternate_allele) = $SEQUENCE=~/[acgtACTG]\(\+([acgtACTG]+)\)/;
			$type='+';
			$ref_seq=~s/\(.+\)//g;
		}
		elsif($SEQUENCE=~/\(\/[acgtACTG]\)/)	# / 
		{
			($reference_allele, $alternate_allele) = $SEQUENCE=~/([acgtACTG])\(\/([acgtACTG])\)/;
			$type='/';
			if(uc($reference_allele) ne uc($alternate_allele))	# / 
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
		# say $SEQUENCE;
		# say "Original sequence:";
		# say $ref_seq;
		# say "Edited sequence:";
		# say $alt_seq;
		
		#         SEQUENCE|**********(+T)****|**********(-T)****|*********A(/T)****
		# reference_allele|                  |            T     |            A     
		# alternate_allele|            T     |                  |            T     
		#          ref_seq|**********    ****|**********  T ****|**********  A ****
		#          alt_seq|**********  T ****|**********    ****|**********  T ****

		my $change=0;
		$change=abs(length($alternate_allele)-length($reference_allele)) if $type eq '+';	
		my $sub_seq=uc(substr($ref_seq,0,$edit_pos+4+length($pam)));
		my $i=rindex($sub_seq,$pam);
		my @locs;
		while( $i>=$len_spacer+1 and $edit_pos-($i-5)<=$CUT_SIZE )
		{
			push @locs,$i;
			$i=rindex($sub_seq,$pam,$i-1);
		}
		
		my $sub_seq_pe3=uc(substr($ref_seq,0,length($ref_seq)-($len_spacer+1)));
		my $i_pe3=index($sub_seq_pe3,rev(complement($pam)));
		my @locs_pe3;
		while( $i_pe3>=0 )
		{
			push @locs_pe3,$i_pe3;
			$i_pe3=index($sub_seq_pe3,rev(complement($pam)),$i_pe3+1);
		}
		
		my @pegrnas;
		# die("There is no PAM and spacer within ${CUT_SIZE}bp from the edit site!!!") unless @locs;
		return @pegrnas unless @locs;
		foreach my $i_pam(@locs)
		{
			my $best_sg=-1;
			my @locs_pe3_filtered;
			foreach my $i_pam_sg(@locs_pe3)
			{
				my $nick_to_peg_dis=($i_pam_sg+3+length($pam))-($i_pam-5);
				next if $nick_to_peg_dis<$MIN_DISGRNA or $nick_to_peg_dis>$MAX_DISGRNA;
				
				push @locs_pe3_filtered,$i_pam_sg;
				my $ind_sg_l=$i_pam_sg+length($pam)+1;
				my $ind_sg_r=$i_pam_sg+length($pam)+$len_spacer;
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
				$best_sg=$locs_pe3_filtered[0];
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
			
			my $pos=$edit_pos-($i_pam-5)+$change;
			$pos-- unless $type eq '/';
			foreach my $len_PBS($MIN_PBS..$MAX_PBS)
			{
				foreach my $len_RT($MIN_RT..$MAX_RT)
				{
					next unless $len_RT>=$pos+$HOMOLOGY;

					my $seq_spacer=substr($ref_seq,$i_pam-1-$len_spacer,$len_spacer);
					my $seq_pam=substr($ref_seq,$i_pam-1,length($pam)+1);
					my $seq_PBS=substr($ref_seq,$i_pam-4-$len_PBS,$len_PBS);
					my $seq_RT=substr($alt_seq,$i_pam-4,$len_RT);
					my $seq_target;
					if($i_pam-1-$len_spacer-$len_up>=0)
					{
						$seq_target=substr($ref_seq,$i_pam-1-$len_spacer-$len_up,$len_target);
					}
					else
					{
						$seq_target=substr($ref_seq,0,$len_target);
					}
					
					$seq_PBS=complement($seq_PBS);
					$seq_PBS=rev($seq_PBS);
					$seq_RT=complement($seq_RT);
					$seq_RT=rev($seq_RT);
					
					$id_peg++;
					push @pegrnas,"$id_peg\t$Strand\t$seq_target\t$seq_spacer\t$seq_pam\t$seq_PBS\t$seq_RT\t".($edit_pos-($i_pam-5))."\t$seq_sg\t$seq_sg_pam\t$nick_to_peg_dis";
					# say OUT "$id_peg\t$seq_target\t$seq_spacer\t$seq_PBS\t$seq_RT\t+".($edit_pos-($i_pam-5));
				}
			}
		}
		return @pegrnas;
		
	}
	else
	{
		die("The input sequence should only contain character of 'acgtACTG+-/()'!!\n");
	}
}


sub complement	
{
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
	join('',reverse(split(//,$_[0])));
}


sub extract_sequence
{
	my $Chromosome;
	my $mid;
	my $len;
	if(@_==3)
	{
		$Chromosome=$_[0];
		$mid=$_[1];
		$len=$_[2];
	}
	else
	{
		$Chromosome='chr1';
		$mid=1014143;
		$len=2;
	}
	my ($start,$end)=($mid-$len,$mid+$len);
	# open(FILE_OUT,">","chr_location.fa");

	# print "$Chromosome:$start:$end\n";

	my $num_colunm=50;
	open(FILE,"<","Data/Reference_Genome/$Chromosome.fa") or die("The input chromosome is not included, and please change to another chromosome!\n");
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
	# print FILE_OUT ">$Chromosome-$start-$end\n$data\n";
	# print ">$Chromosome\t$start\t$end\n$data\n";
	close FILE;
	return $data;
}