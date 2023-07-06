#!/usr/bin/env perl
use strict;
# use warnings;
use 5.010;


# Example for substitution: CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTGGCGCGA(/T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC
# Example for insertion: CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTG(+ATT)GCGCGAGGCCGCCTGGCAACTCTGCGACTACTACCTGCC
# Example for deletion: CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGC(-GCTGGCGCGA)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC

#  xuhao =====================
# my $SEQUENCE='aaaaCGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTGGCGCGA(+T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC';
# # my $SEQUENCE='aaaaCGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTGGCGCGA(-T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC';
# # my $SEQUENCE='aaaaCGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTGGCGCGA(/T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC';
# my $PAM='NGG';	#PAM
# my $CUT_SIZE=10;	#maximum edit-to-nick distance 
# my ($MIN_DISGRNA,$MAX_DISGRNA)=(0,100);	#minimum/maximum nick-to-nick distance
# my ($MAX_PBS,$MAX_RT)=(17,24);
# my ($MIN_PBS,$MIN_RT)=(8,8);
# my $HOMOLOGY=5;	#minimum downstream HOMOLOGY
# my $time=1234567;
# print "$SEQUENCE\n";
# chdir '../';
my $SEQUENCE=$ARGV[0];
my $PAM=$ARGV[1];
my $CUT_SIZE=$ARGV[2];
my $MIN_DISGRNA=$ARGV[3];
my $MAX_DISGRNA=$ARGV[4];
my $MAX_PBS=$ARGV[5];
my $MAX_RT=$ARGV[6];
my $MIN_PBS=$ARGV[7];
my $MIN_RT=$ARGV[8];
my $HOMOLOGY=$ARGV[9];
my $time=$ARGV[10];
#  xuhao =====================



##
my $type='User';
my $pam=substr($PAM,1);	# remove N, for convenience to match
my $len_spacer=20;	#spacer  20bp
my $len_target=47;	#Total 47 bps = 4 bp neighboring sequence + 20 bp protospacer + 3 bp NGG PAM+ 20 bp neighboring sequence
my $len_up=4;	#NBT spacer upper
my $dir_out="Temp";
mkdir $dir_out unless -e $dir_out;
open(OUT,">","$dir_out/Sequence.request.$type.$time.txt") or die($!);
say OUT "PegRNAID\tStrand\tTarget\tSpacer\tPAM\tPBS\tRT\tEditToNickDistance\tsgRNASpacer\tsgRNAPAM\tNickToNickDistance";

my $user_seq_minus='';		# negtive strand
my ($up, $mid, $down)=("","","");
die("The input sequence should only contain character of 'acgtACTG+-/()'!\n") unless $SEQUENCE=~/^[acgtACTG+-\/()]+$/;
if($SEQUENCE=~/^[acgtACTG]*\(-[acgtACTG]+\)[acgtACTG]*$/)	#  -  deletion
{
	($up, $mid, $down) = $SEQUENCE=~/^([acgtACTG]*)\(-([acgtACTG]+)\)([acgtACTG]*)$/;
	$user_seq_minus=rev(complement($down)).'(-'.rev(complement($mid)).')'.rev(complement($up));
}
elsif($SEQUENCE=~/^[acgtACTG]*\(\+[acgtACTG]+\)[acgtACTG]*$/)	#  +  insertion
{
	($up, $mid, $down) = $SEQUENCE=~/^([acgtACTG]*)\(\+([acgtACTG]+)\)([acgtACTG]*)$/;
	$user_seq_minus=rev(complement($down)).'(+'.rev(complement($mid)).')'.rev(complement($up));
}
elsif($SEQUENCE=~/^[acgtACTG]+\(\/[acgtACTG]\)[acgtACTG]*$/)	#  /  substitution
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
	Example for substitution: CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTGGCGCGA(/T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC
	Example for insertion: CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTG(+ATT)GCGCGAGGCCGCCTGGCAACTCTGCGACTACTACCTGCC
	Example for deletion: CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGC(-GCTGGCGCGA)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC
	';
	die("The format of the input sequence with the intended edit is incorrect!! Please check and input a correct one!!\n");
}
# die('There should be at least one nucleotide before the parentheses') if($up eq "");
# die('There should be at least one nucleotide after the parentheses') if($down eq "");

my $count_variants=0;
my $count_rnas=0;

my @pegrnas=search_candicante($SEQUENCE,0,'+');	# positive strand
if(@pegrnas)
{
	$count_rnas+=@pegrnas;
	foreach my $rna(@pegrnas)
	{
		say OUT "$rna";
	}
}

my $id=0;
$id=(split(/\t/,$pegrnas[-1]))[0] if @pegrnas;	#get the last pegRNA ID in + strand
my @pegrnas_minus=search_candicante($user_seq_minus,$id,'-');	# - strand
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
	my $id_peg=$_[1];	# start 0
	my $Strand=$_[2];	# strand + / -
	my $ref_seq=$SEQUENCE;
	if($SEQUENCE=~/^[acgtACTG+-\/()]+$/)
	{
		my $edit_pos =index($SEQUENCE,'(')-1;	# 0 for the first bp
#		die("The current length of upstream sequence is ".($edit_pos+1)."bp and should be at least ".($len_spacer-3)."bp.") if $edit_pos+1<=$len_spacer-3;	# must >= 17bp upper stand
		
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
		# print "$change\n";
		my $sub_seq=uc(substr($ref_seq,0,$edit_pos+4+length($pam)));
		my $i=rindex($sub_seq,$pam);	# locate 2nd base G of PAM
		my @locs;	#pegRNA location 2nd base G of PAM
		while( $i>=$len_spacer+1 and $edit_pos-($i-5)<=$CUT_SIZE )	#$i>=$len_spacer+1 PAM 20bp spacer；$i-5 - index of cut，$edit_pos - edition position
		{
			push @locs,$i;
			$i=rindex($sub_seq,$pam,$i-1);		
		}
		# print "@locs \n";
		
		my $sub_seq_pe3=uc(substr($ref_seq,0,length($ref_seq)-($len_spacer+1)));
		my $i_pe3=index($sub_seq_pe3,rev(complement($pam)));	# location of first C of CCN / CN
		my @locs_pe3;	#PE3 sgRNA location of first C of CCN / CN
		while( $i_pe3>=0 )
		{
			push @locs_pe3,$i_pe3;
			$i_pe3=index($sub_seq_pe3,rev(complement($pam)),$i_pe3+1);
		}
		
		my @pegrnas;
		# die("There is no PAM and spacer within ${CUT_SIZE}bp from the edit site!!!") unless @locs;	# no PAM
		return @pegrnas unless @locs;
		foreach my $i_pam(@locs)
		{
			my $best_sg=-1;	#
			my @locs_pe3_filtered;	# ok PE3 sgRNA location of first C of CCN / CN
			foreach my $i_pam_sg(@locs_pe3)
			{
				my $nick_to_peg_dis=($i_pam_sg+3+length($pam))-($i_pam-5);
				next if $nick_to_peg_dis<$MIN_DISGRNA or $nick_to_peg_dis>$MAX_DISGRNA;
				
				push @locs_pe3_filtered,$i_pam_sg;
				my $ind_sg_l=$i_pam_sg+length($pam)+1;	#sgRNA 20nt spacer the first from left
				my $ind_sg_r=$i_pam_sg+length($pam)+$len_spacer;	#sgRNA 20nt spacer the first from right
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
				$best_sg=$locs_pe3_filtered[0];	# primedesign, PE3 annotation at a distance as close to 75 bp away
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
			
			my $pos=$edit_pos-($i_pam-5)+$change;	#pos position of the first base from right,$edit_pos-($i_pam-5) is the index of edit point relative to cut point, +1, +2, ...
			$pos-- unless $type eq '/';
			foreach my $len_PBS($MIN_PBS..$MAX_PBS)
			{
				foreach my $len_RT($MIN_RT..$MAX_RT)
				{
					next unless $len_RT>=$pos+$HOMOLOGY;	#RT at least 5bp betweeb edit point to the first from the right

					my $seq_spacer=substr($ref_seq,$i_pam-1-$len_spacer,$len_spacer);
					my $seq_pam=substr($ref_seq,$i_pam-1,length($pam)+1);
					my $seq_PBS=substr($ref_seq,$i_pam-4-$len_PBS,$len_PBS);
					my $seq_RT=substr($alt_seq,$i_pam-4,$len_RT);
					my $seq_target;
					if($i_pam-1-$len_spacer-$len_up>=0)
					{
						$seq_target=substr($ref_seq,$i_pam-1-$len_spacer-$len_up,$len_target);	# seq_target: Spacer upper+Spacer+PAM+Spacer downstream. Total 47 bps = 4 bp neighboring sequence + 20 bp protospacer + 3 bp NGG PAM+ 20 bp neighboring sequence
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
