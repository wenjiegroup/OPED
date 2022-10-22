#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

#更新前字段
#my $alleleID=$ARGV[0];
#my $pam_type=$ARGV[1];
#my $direction=$ARGV[2];


#一个简单的查表操作（文件在Database目录中,final_pegRNA开头的两个文件），用户输入AlleleID、 选择PAM类型（NGG或NG）、选择编辑方向（install 或 correct），然后脚本会打开对应文件，找到对应AlleleID，返回模型设计的最优的top10个pegRNA

##以下几个参数为用户输入的，否则采用下面的默认参数
#my $query_type='AlleleID';	# AlleleID	GeneID	GeneSymbol	HGNC_ID 共支持这4种查询类型， 需要选择一种
#my $query_item='929884';	# 在确定$query_type查询类型的具体查询的项，这里让用户输入
#my $pam_type='NGG';	#NGG或NG
#my $direction='install';	#install or correct pathogenic human genetic variants from the ClinVar database

my $query_type=$ARGV[0];
my $query_item=$ARGV[1];
my $pam_type=$ARGV[2];
my $direction=$ARGV[3];

##此部分为固定的参数
my $top=10; 
#my $dir_in="Database";
#via hjs 四月
my $dir_in="Data";


open(IN,"<","$dir_in/final_pegRNA.$direction.top$top.txt") or die("The input parameter is wrong and the file does not exit!\n");
my $header=<IN>;
$header=~s/[\r\n]$//;
my @temp=split(/\t/, $header);
my $index=-1;
foreach my $i(0..$#temp)
{
	if($temp[$i] eq $query_type)
	{
		$index=$i;
		last;
	}
}
die("The input type of query is incorrect. \nOnly support query for AlleleID or GeneID or GeneSymbol or HGNC_ID\n") if $index<0;

my @pegRNAs;
while(<IN>)
{
	s/[\r\n]+$//;
	my @temp=split(/\t/);
	if($temp[$index] eq $query_item and $temp[-1] eq $pam_type)
	{
		push @pegRNAs, "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\t$temp[11]\t$temp[13]\t$temp[15]\t$temp[16]\t$temp[17]\t$temp[18]\t$temp[20]";
	}
}
if(@pegRNAs)
{
	open(OUT,">","$dir_in/output_of_search.txt") or die($!);
	my @temp=split(/\t/,$header);
	say OUT "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\t$temp[11]\t$temp[13]\t$temp[15]\t$temp[16]\t$temp[17]\t$temp[18]\t$temp[20]";
	# say "Select top $top but get ".scalar(@pegRNAs) unless scalar(@pegRNAs)==$top;
	foreach (@pegRNAs)
	{
		say OUT;
	}
}
else
{
	die("The optimal pegRNA(PAM=$pam_type) to $direction pathogenic human genetic variants of query $query_type:$query_item does not exit in the OPEDVar database!\n");
}



