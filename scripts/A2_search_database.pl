#!/usr/bin/env perl
use DBI;
use strict;
use warnings;
use 5.010;


#  xuhao =====================
# my $queryType='AlleleID';	# AlleleID	GeneID	GeneSymbol	HGNC_ID 
# my $query_item='9298d84';
# my $PAM='NGG';	#NGG / NG
# my $direction='install';	#install or correct pathogenic human genetic variants from the ClinVar database
# my $time='20230405080000BJT';
# chdir '../';
my $queryType=$ARGV[0];
my $query_item=$ARGV[1];
my $PAM=$ARGV[2];
my $direction=$ARGV[3];
my $time=$ARGV[4];
#  xuhao =====================

my $top=10; 
#my $dir_in="Database";
my $dir_in="Data";
my $dir_out="Temp";
mkdir $dir_out unless -e $dir_out;

my $driver   = "SQLite"; 
my $database = "Data/OPEDVar.db";
my $dsn = "DBI:$driver:dbname=$database";
my $dbh = DBI->connect($dsn, "", "", { RaiseError => 1 }) or die $DBI::errstr;
print "Opened database successfully\n";


my $stmt = qq(SELECT AlleleID, Type, GeneID, GeneSymbol, HGNC_ID, Chromosome, Start, Stop, ReferenceAllele, AlternateAllele,
	Strand, Spacer, PBS, RT, EditToNickDistance, sgRNASpacer, NickToNickDistance, PAM
	FROM `final_pegRNA.$direction.top10`
	WHERE $queryType = '$query_item' AND PAM = '$PAM';);
my $sth = $dbh->prepare( $stmt );
my $rv = $sth->execute() or die $DBI::errstr;
if($rv < 0){
	print $DBI::errstr;
}

open(OUT,">","$dir_out/OPEDVar.request.$time.txt") or die($!);
say OUT "AlleleID\tType\tGeneID\tGeneSymbol\tHGNC_ID\tChromosome\tStart\tStop\tReferenceAllele\tAlternateAllele\tStrand\tSpacer\tPBS\tRT\tEditToNickDistance\tsgRNASpacer\tNickToNickDistance\tPAM";

my $n = 0;

while(my @row = $sth->fetchrow_array()) {
	$n = $n + 1;
	say OUT "$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]\t$row[8]\t$row[9]\t$row[10]\t$row[11]\t$row[12]\t$row[13]\t$row[14]\t$row[15]\t$row[16]\t$row[17]"
}

if ($n eq 0){
	die("NO RECORDS FOUND: The optimized pegRNA to $direction pathogenic human genetic variants of query does not exit in the OPEDVar database!\n* $queryType=$query_item\n* (PAM=$PAM)\n");
}

print "Operation done successfully\n";
$dbh->disconnect();




