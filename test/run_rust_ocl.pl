use strict;
use warnings;
use Cwd qw(getcwd);


my $currentdir =  getcwd();

my $rustdir = "../rust/sa_opencl";
chdir $rustdir;
system("cargo build --release");

chdir $currentdir;

my $seqdir  ="inputs/";
opendir(DIR,$seqdir);
my @allfiles = grep(/\.1\.fas$/,readdir(DIR));
closedir(DIR);
my $outdir = "res_rust_ocl/";

if(!-d $outdir){
	mkdir($outdir) or die;
}

#my @options=("local","glocal","global");
my @options=("local","glocal","global");
foreach my $oo(@options){
	foreach my $ss(@allfiles){
		my $infile1 = $seqdir."/".$ss;
		my $infile2 = $seqdir."/".$ss;
		$infile2 =~ s/\.1\.fas$/.2.fas/;
		my $outname = $outdir."/".$ss.".res.$oo";
		my @res = `$rustdir/target/release/sa_opencl -$oo $infile1 $infile2`;
		open(OUT,"> $outname");
		foreach my $ll(@res){
			print OUT $ll;
		}
		close(OUT);
	}
}