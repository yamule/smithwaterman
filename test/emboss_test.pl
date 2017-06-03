use strict;
use warnings;


# kiga muitara kireini suru wa,,,

my $embosswater = "EMBOSS-6.6.0/emboss/water ";#path to the emboss water http://emboss.sourceforge.net/
my $sourcefas = "hsa_sprot.fasta";# multiple fasta file 

my ($heads,$aas) = getFasta($sourcefas);
my $recheck = 0;
for(my $i = 0;$i < 100;$i++){
	my $sa = int(rand()*($#{$heads}+1));
	my $sb = int(rand()*($#{$heads}+1));
	
	my $fasa = ${$aas}[$sa];
	my $fasb = ${$aas}[$sb];
	
	if($recheck == 0){
		open(OUT,"> s1.fas");
		print OUT ">s1\n";
		print OUT $fasa;
		close(OUT);
		open(OUT,"> s2.fas");
		print OUT ">s2\n";
		print OUT $fasb;
		close(OUT);
	}
	
	system("$embosswater -asequence s1.fas -bsequence s2.fas -gapopen 10.0 -gapextend 0.5  -outfile testout.dat ");
	open(IN, "testout.dat");
	my $resa =  "";
	my $resb =  "";
	#get aligned sequences from emboss water 
	while(my $ss = <IN>){
		if($ss =~ /^s1[\s]+[0-9]+[\s]+([^\s]+)/){
			$resa .= $1;
		}
		if($ss =~ /^s2[\s]+[0-9]+[\s]+([^\s]+)/){
			$resb .= $1;
		}
	}
	close(IN);
	
	
	# scripts must print 
	# #score
	# <name of the first sequence> <alignedseq>
	# <name of the second sequence> <alignedseq>
	# in print stdout 
	#system("python3 smithwaterman.py s1.fas s2.fas > testout.dat ");
	system("perl smithwaterman_run.pl s1.fas s2.fas > testout.dat ");
	open(IN, "testout.dat");
	my $ssa =  "";
	my $ssb =  "";
	my $score = 0;
	while(my $ss = <IN>){
		if($ss =~ /^#score[^0-9]+([0-9]+)/){
			$score = $1;
		}
		if($ss =~ /^s1[\s]+([^\s]+)/){
			$ssa .= $1;
		}
		if($ss =~ /^s2[\s]+([^\s]+)/){
			$ssb .= $1;
		}
	}
	close(IN);
	
	
	
	my $ca = $ssa;
	my $cb = $ssb;
	$ca =~ s/[^A-Z]//g;
	$cb =~ s/[^A-Z]//g;
	$fasa =~ s/[^A-Z]//g;
	$fasb =~ s/[^A-Z]//g;
	if($ca eq $fasa && $cb eq $fasb){
		print "full length OK\n";
	}else{
	
		print "//--------------------------- full length ng\n";
		
		print ">s1\n";
		print $ca."\n";
		print ">s2\n";
		print $cb."\n";
		print ">s1\n";
		print $fasa."\n";
		print ">s2\n";
		print $fasb."\n";
		#last;
	}
	
	
	
	
	my @aaa = split(//,$ssa);
	my @bbb = split(//,$ssb);
	
	my $st = -1;
	my $en = -2;
	for(my $aa = 0;$aa <= $#aaa;$aa++){
		if($aaa[$aa] =~ /[A-Z]/ && $bbb[$aa] =~ /[A-Z]/){
			$st = $aa;
			last;
		}
	}
	for(my $aa = $#aaa;$aa >= 0;$aa--){
		if($aaa[$aa] =~ /[A-Z]/ && $bbb[$aa] =~ /[A-Z]/){
			$en = $aa;
			last;
		}
	}
	
	$ssa = join("",@aaa[$st..$en]);
	$ssb = join("",@bbb[$st..$en]);
	if($ssa eq $resa && $ssb eq $resb){
		print $score."\n";
		print "alignment with water OK\n";
	}else{
		print "emboss: ".$resa."\n";
		print "script: ".$ssa."\n";
		print "\n";
		print "emboss: ".$resb."\n";
		print "script: ".$ssb."\n";
		if(1 == 0){
			print "//---------------------------\n";
			
			print ">s1\n";
			print $fasa."\n";
			print ">s2\n";
			print $fasb."\n";
		}
		last;
	}
	if($recheck == 1){
		last;
	}
}




sub getFasta{
	my $filename = $_[0];
	
	my @head;
	my @seq;
	open(IN,$filename);
	while(my $ss = <IN>){
		if($ss =~ />/){
			push(@head,$ss);
			push(@seq,"");
		}else{
			$seq[$#head] .= $ss;
		}
		
		
	}
	close(IN);
	return \@head,\@seq;
}
