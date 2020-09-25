package sw;

use strict;
use warnings;


# public domain
# no warranty
# author: yamule (https://github.com/yamule)
# usage:
# require "smithwaterman.pl";
# my ($res1,$res2,$score) = sw::align("QUERYSEQ", "TEMPLATESEQ");
# print ">seq1\n".join("",@{$res1})."\n";
# print ">seq2\n".join("",@{$res2})."\n";
# Output:
# >seq1
# --------QXERYSEQ
# >seq2
# TEMPLATE-----SEQ


my @sw_mat=(
#https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by sw_matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *",
"A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 ",
"R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 ",
"N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 ",
"D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 ",
"C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 ",
"Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 ",
"E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ",
"G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 ",
"H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 ",
"I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 ",
"L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 ",
"K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 ",
"M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 ",
"F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 ",
"P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 ",
"S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 ",
"T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 ",
"W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 ",
"Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 ",
"V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 ",
"B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 ",
"Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ",
"X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 ",
"* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 "
);

my %subscore;
my @hcode;
for(my $ii = 0;$ii <= $#sw_mat;$ii++){
	$sw_mat[$ii] =~ s/^[\s]+//g;
	$sw_mat[$ii] =~ s/[\s]+$//g;
	my @ar = split(/[\s]+/,$sw_mat[$ii]);
	if($ii == 0){
		@hcode = @ar;
	}else{
		my $dcode = $ar[0];
		my %tmp;
		$subscore{$dcode} = \%tmp;
		for(my $jj = 1;$jj <= $#ar;$jj++){
			${$subscore{$dcode}}{$hcode[$jj-1]} = $ar[$jj];
		}
	}
}




#my $seq1 = "RRHRSEDCGG--------GPRS--LSRGLPCKKAATEGS---------SEKTVLDSKPSVPTTSEGGPELELQIPELPLDSNEFWVHEGCILWANGIYLVCGRLYGLQEALEIAREMKCSHCQEAGATLGCYN--KGCSFRYHYPCAIDADCLLHEENFSVRCPKHKPPLPCPLPPLQNKTAKGSLSTEQSERG--";
#my $seq2 = "----------MKLAFLFLGPMALLLLAGYGCVLGASSGNLRTFVGCAVREFTFLAKKPG----CRG----------LRITTDACWGR--CETWEKPI---------LEPPYIEAHHRVCTYNETKQVTVKLPNCAPGVDPFYTYPVAIRCDC-------------------------------GACSTATTECETI";
my $seq1 = "MSIGVPIKVLHEAEG-----HIVTCETNTGEVYRGKLIEAEDNMNCQMSNI-----TVTYRDGRVAQLEQVYIRGSK-----------IRFLILPDMLKNAPMLKSMKNKNQGSGAGRGKAAILKAQVAARGRGRGMGRGNIFQKRR";
my $seq2 = "MAVAWGIGFLHSVSQLAFAVHLPFCGPN--EV---------DSFYCDLPRVIKLACTDTYR------LDIMVIANSGVLTVCSFVLLIISYTIILMTIQHCPLDKSSK-----------ALSTLTAHIT------------------";

my @aa1 = split(//,$seq1);
my @aa2 = split(//,$seq2);
my $score = 0;
my $flag = 0;
for(my $ii = 0;$ii <= $#aa1;$ii++){
	
	if($aa1[$ii] !~ /-/ && $aa2[$ii] !~ /-/){
		$flag = 1;
		$score+=${$subscore{$aa1[$ii]}}{$aa2[$ii]};
	}else{
		if($flag == 1){
			if($aa1[$ii-1] =~ /-/ || $aa2[$ii-1] =~ /-/){
				$score -= 0.5;
			}else{
				$score -= 10.0;
			}
		}
	}
}

print $score."\n";
