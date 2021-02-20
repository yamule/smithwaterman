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
my $outdir = "res_rust_ocl_list/";

if(!-d $outdir){
	mkdir($outdir) or die;
}

#my @options=("local","glocal","global");
my @options=("local","glocal","global");
foreach my $oo(@options){
	my $outname = $outdir.$oo."tmp.dat";
	my @res = `$rustdir/target/release/sa_opencl -$oo -list file_list.txt `;
	open(OUT,"> $outname");
	foreach my $ll(@res){
		print OUT $ll;
	}
	close(OUT);
	open(IN,"file_list.txt");
	my @file_basename;
	while(my $ss = <IN>){
		if($ss =~ /\/([^\t]+)\t/){
			push(@file_basename,$1);
		}
	}
	close(IN);
	my $buff = "";
	open(IN,$outname);
	while(my $ss = <IN>){
		if($ss =~ /#score/){
			if(length($buff) > 0){
				my $bname = shift(@file_basename);
				open(OOUT,"> $outdir$bname.res.$oo");
				print OOUT $buff;
				close(OOUT);
			}
			$buff = "";
		}
		$buff .= $ss;
	}
	close(IN);
	
	if(length($buff) > 0){
		my $bname = shift(@file_basename);
		open(OOUT,"> $outdir$bname.res.$oo");
		print OOUT $buff;
		close(OOUT);
	}
}