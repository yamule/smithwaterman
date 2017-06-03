package sw;

use strict;
use warnings;


# public domain
# no warranty
# author: yamule (https://github.com/yamule)
# usage:
# my ($res1,$res2) = sw::align("QUERYSEQ", "TEMPLATESEQ");
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







sub align{
	my $seq1 = uc $_[0];
	my $seq2 = uc $_[1];
	
	my @reta;
	my @retb;
	my $MAT_ = 0;
	my $GAPA_ = 1;
	my $GAPB_ = 2;

	my $go = -10;
	my $ge = -0.5;
	
	
	$seq1 =~ s/[^A-Za-z]//g;
	$seq2 =~ s/[^A-Za-z]//g;
	
	$seq1 =~ s/[BJOUXZa-z]/X/g;
	$seq2 =~ s/[BJOUXZa-z]/X/g;
	
	my @aa1 = split(//,$seq1);
	my @aa2 = split(//,$seq2);
	
	
	my @dpsw_mat;
	for(my $a = 0;$a <= $#aa1+1;$a++){
		my @tmp;
		push(@dpsw_mat,\@tmp);
		for(my $b = 0;$b <= $#aa2+1;$b++){
			my @tmp2;
			push(@tmp,\@tmp2);
		}
	}
	
	
	for(my $a = 0;$a <= $#aa1+1;$a++){
		for(my $b = 0;$b <= $#aa2+1;$b++){
			if($a == 0 || $b == 0){
				${${$dpsw_mat[$a]}[$b]}[$MAT_] = 0;
				${${$dpsw_mat[$a]}[$b]}[$GAPA_] = 0;
				${${$dpsw_mat[$a]}[$b]}[$GAPB_] = 0;
			}
		}
	}
	for(my $a = 1;$a <= $#aa1+1;$a++){
		for(my $b = 1;$b <= $#aa2+1;$b++){
			${${$dpsw_mat[$a]}[$b]}[$MAT_] = getMax(getMax3(${${$dpsw_mat[$a-1]}[$b-1]}[$MAT_]
			,${${$dpsw_mat[$a-1]}[$b-1]}[$GAPA_]
			,${${$dpsw_mat[$a-1]}[$b-1]}[$GAPB_]) + ${$subscore{$aa1[$a-1]}}{$aa2[$b-1]},0);
			
			
			
			${${$dpsw_mat[$a]}[$b]}[$GAPA_] = getMax3(${${$dpsw_mat[$a]}[$b-1]}[$MAT_]+$go,${${$dpsw_mat[$a]}[$b-1]}[$GAPA_]+$ge,0);
			${${$dpsw_mat[$a]}[$b]}[$GAPB_] =  getMax3(${${$dpsw_mat[$a-1]}[$b]}[$MAT_]+$go,${${$dpsw_mat[$a-1]}[$b]}[$GAPB_]+$ge,0);
		}
	}
	
	
	
	my $starta = -1;
	my $startb = -1;
	my $startcode = -1;
	my $maxscore = -1;
	for(my $a = 1;$a <= $#aa1+1;$a++){
		for(my $b = 1;$b <= $#aa2+1;$b++){
			for(my $ii = 0;$ii < 3;$ii++){
				if(${${$dpsw_mat[$a]}[$b]}[$ii] > $maxscore){
					$maxscore = ${${$dpsw_mat[$a]}[$b]}[$ii];
					$starta = $a;
					$startb = $b;
					$startcode = $ii;
				}
			}
		}
	}
	
	my $currentcode = $startcode;
	if($currentcode == $MAT_){#MATˆÈŠO‚È‚¢‚Í‚¸‚¾‚ªB
		my $sa = $starta;
		my $sb = $startb;
		
		for(my $ii = $#aa1;$ii >= $sa;$ii--){
			push(@reta,$aa1[$ii]);
			push(@retb,"-");
			
		}
		for(my $ii = $#aa2;$ii >= $sb;$ii--){
			push(@retb,$aa2[$ii]);
			push(@reta,"-");
		}
		
		
		push(@reta,$aa1[$sa-1]);
		push(@retb,$aa2[$sb-1]);
	}else{
		die;
	}
	while(1==1){
		my $ms = -1;
			
		if($currentcode == $MAT_){
			my $sa = $starta-1;
			my $sb = $startb-1;
			
			my $cd = -1;
			for(my $ii = 0;$ii < 3;$ii++){
				if(${${$dpsw_mat[$sa]}[$sb]}[$ii] > $ms){
					$ms = ${${$dpsw_mat[$sa]}[$sb]}[$ii];
					$cd = $ii;
				}
			}
			
			$currentcode = $cd;
			$starta = $sa;
			$startb = $sb;
			
		}elsif($currentcode == $GAPA_){
			
			my $sa = $starta;
			my $sb = $startb-1;
			my $cd = -1;
			
			$ms = ${${$dpsw_mat[$sa]}[$sb]}[$MAT_];
			$cd = $MAT_;
			
			
			if(${${$dpsw_mat[$sa]}[$sb]}[$GAPA_]+$ge > $ms+$go){
				$ms = ${${$dpsw_mat[$sa]}[$sb]}[$GAPA_];
				$cd = $GAPA_;
			}
			
			
			$currentcode = $cd;
			$starta = $sa;
			$startb = $sb;
		}elsif($currentcode == $GAPB_){
			
			my $sa = $starta-1;
			my $sb = $startb;
			my $cd = -1;
			
			
			$ms = ${${$dpsw_mat[$sa]}[$sb]}[$MAT_];
			$cd = $MAT_;
			if(${${$dpsw_mat[$sa]}[$sb]}[$GAPB_]+$ge > $ms+$go){
				$ms = ${${$dpsw_mat[$sa]}[$sb]}[$GAPB_];
				$cd = $GAPB_;
			}
			
			$currentcode = $cd;
			$starta = $sa;
			$startb = $sb;
		}
		if($ms <= 0){
			if($currentcode == $MAT_){
				$starta++;
				$startb++;
			}elsif($currentcode == $GAPA_){
				$startb++;
			}elsif($currentcode == $GAPB_){
				$starta++;
			}
			last;
		}
		if($currentcode == $MAT_){
			push(@reta,$aa1[$starta-1]);
			push(@retb,$aa2[$startb-1]);
		}elsif($currentcode == $GAPA_){
			push(@reta,"-");
			push(@retb,$aa2[$startb-1]);
		}elsif($currentcode == $GAPB_){
			push(@reta,$aa1[$starta-1]);
			push(@retb,"-");
		}
		
	}
	
	for(my $ii = $starta-2;$ii >= 0;$ii--){
		push(@reta,$aa1[$ii]);
		push(@retb,"-");
		
	}
	for(my $ii = $startb-2;$ii >= 0;$ii--){
		push(@retb,$aa2[$ii]);
		push(@reta,"-");
		
	}
	
	
	
	my @retaa = reverse(@reta);
	my @retbb = reverse(@retb);
	return \@retaa,\@retbb,$maxscore;
}



sub getMax{
	my $a = $_[0];
	my $b = $_[1];
	if($a >= $b){
		return $a;
	}
	return $b;
	
}





sub getMax3{
	my $a = $_[0];
	my $b = $_[1];
	my $c = $_[2];
	if($a >= $b){
		if($a >= $c){
			return $a;
		}else{
			return $c;
		}
	}
	if($b >= $c){
		return $b;
	}
	return $c;
	
}



1;