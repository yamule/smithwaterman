use strict;
use warnings;

my $classdir = "javaclass/";
mkdir($classdir);
system("javac ../SmithWaterman.java -d $classdir");



my $seqdir  ="inputs/";
opendir(DIR,$seqdir);
my @allfiles = grep(/\.1\.fas$/,readdir(DIR));
closedir(DIR);
my $outdir = "res_java/";
if(!-d $outdir){
	mkdir($outdir) or die;
}

foreach my $ss(@allfiles){
	my $infile1 = $seqdir."/".$ss;
	my $infile2 = $seqdir."/".$ss;
	$infile2 =~ s/\.1\.fas$/.2.fas/;
	my $outname = $outdir."/".$ss.".res";
	my @res = `java -cp $classdir SmithWaterman $infile1 $infile2`;
	open(OUT,"> $outname");
	foreach my $ll(@res){
		print OUT $ll;
	}
	close(OUT);
}

