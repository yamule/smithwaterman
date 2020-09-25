use strict;
use warnings;


my $embossneedle = "needle";#path to the emboss water http://emboss.sourceforge.net/
my $embosswater = "water";#path to the emboss water http://emboss.sourceforge.net/
my $sourcefas = "../../hsa_sprot.fas";# multiple fasta file 

my $num_trials = 100; #how many alignments will be made.

my ($heads,$aas) = getFasta($sourcefas);
my $indir = "inputs";
my $outdir = "emboss_results";
mkdir($indir);
mkdir($outdir);
for(my $i = 0;$i < $num_trials;$i++){
	my $sa = int(rand()*($#{$heads}+1));
	my $sb = int(rand()*($#{$heads}+1));
	
	my $fasa = ${$aas}[$sa];
	my $fasb = ${$aas}[$sb];
	
	
	my $infile1 = $indir."/seq".$i.".1.fas";
	my $infile2 = $indir."/seq".$i.".2.fas";
	my $outfilename = $outdir."/needle_res".$i.".dat";
	my $outfilename2 = $outdir."/needle_glocal_res".$i.".dat";
	my $outfilename3 = $outdir."/res".$i.".dat";
	
	
	open(OUT,"> $infile1");
	print OUT ">s1\n";
	print OUT $fasa;
	close(OUT);
	open(OUT,"> $infile2");
	print OUT ">s2\n";
	print OUT $fasb;
	close(OUT);
	
	system("$embossneedle -endweight Y  -asequence $infile1 -bsequence $infile2 -gapopen 10.0 -gapextend 0.5 -outfile $outfilename -datafile EBLOSUM62 -sprotein1 -sprotein2 ");
	system("$embossneedle  -asequence $infile1 -bsequence $infile2 -gapopen 10.0 -gapextend 0.5 -outfile $outfilename2  -datafile EBLOSUM62  -sprotein1 -sprotein2 ");
	system("$embosswater -asequence $infile1 -bsequence $infile2 -gapopen 10.0 -gapextend 0.5  -outfile $outfilename3  -datafile EBLOSUM62  -sprotein1 -sprotein2 ");
	

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
