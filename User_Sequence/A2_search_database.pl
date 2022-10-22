#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;


#一个简单的查表操作（文件在Database目录中,final_pegRNA开头的两个文件），用户输入AlleleID、 选择PAM类型（NGG或NG）、选择编辑方向（install 或 correct），然后脚本会打开对应文件，找到对应AlleleID，返回模型设计的最优的top10个pegRNA

##以下几个参数为用户输入的，否则采用下面的默认参数
my $alleleID=$ARGV[0]
my $pam_type=$ARGV[1]
my $direction=$ARGV[2]

#my $alleleID='929884';	# AlleleID of pathogenic human genetic variants from the ClinVar database
#my $pam_type='NGG';	#NGG或NG
#my $direction='install';	#install or correct pathogenic human genetic variants from the ClinVar database

##此部分为固定的参数
my $top=10; 
#my $dir_in="Database";
#myself
my $dir_in="User_Sequence/Database"

open(IN,"<","$dir_in/final_pegRNA.$direction.top$top.txt") or die("$!The input parameter is wrong and the file does not exit!");
my $header=<IN>;
$header=~s/[\r\n]$//;
my @pegRNAs;
while(<IN>)
{
	s/[\r\n]+$//;
	my @temp=split(/\t/);
	if($temp[0] eq $alleleID and $temp[-1] eq $pam_type)
	{
		push @pegRNAs, "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[10]\t$temp[11]\t$temp[12]\t$temp[-1]";
	}
}
if(@pegRNAs)
{
	open(OUT,">","$dir_in/output_of_search.txt") or die($!);
	my @temp=split(/\t/,$header);
	say OUT "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[10]\t$temp[11]\t$temp[12]\t$temp[-1]";
	say "Select top $top but get ".scalar(@pegRNAs) unless scalar(@pegRNAs)==$top;
	foreach (@pegRNAs)
	{
		say OUT;
	}
}
else
{
	die("The pegRNA(PAM=$pam_type) to $direction pathogenic human genetic variants AlleleID:$alleleID from the ClinVar database does not exit!");
}



