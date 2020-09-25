use strict;
use warnings;


require "../smithwaterman.pl";
#make formatted output of smithwaterman.pl

my $seqdir  ="inputs/";
opendir(DIR,$seqdir);
my @allfiles = grep(/\.1\.fas$/,readdir(DIR));
closedir(DIR);
my $outdir = "res_pl";
if(!-d $outdir){
	mkdir($outdir) or die;
}

foreach my $ss(@allfiles){
	my $infile1 = $seqdir."/".$ss;
	my $infile2 = $seqdir."/".$ss;
	$infile2 =~ s/\.1\.fas$/.2.fas/;
	my $outname = $outdir."/".$ss.".res";
	alignFile($infile1,$infile2,$outname);
}

sub alignFile{
	my $file1 = $_[0];
	my $file2 = $_[1];
	my $outfile = $_[2];
	if(-f $outfile){
		die;
	}
	my %s1 = %{${loadFasta($file1)}[0]};
	my %s2 = %{${loadFasta($file2)}[0]};

	my ($res1,$res2,$score) = sw::align($s1{"seq"},$s2{"seq"});
	open(OUT,">".$outfile);
	print OUT "#score: $score\n";
	print OUT $s1{"name"}." ".join("",@{$res1})."\n";
	print OUT $s2{"name"}." ".join("",@{$res2})."\n";
	close(OUT);
}



sub loadFasta{
	my $filename = $_[0];
	my @ret;
	my %t;
	my $currenthash = \%t;
	$t{"seq"} = "";
	$t{"name"} = "";
	$t{"desc"} = "";
	open(IN,$filename);
	while(my $ss = <IN>){
		if($ss =~ />/){
			my $name = "";
			my $desc = "";
			my %tm;
			push(@ret,\%tm);
			if($ss =~ /^[\s]*>[\s]*([^\s]+)/){
				$name = $1;
			}
			if($ss =~ /^[\s]*>[\s]*([^\s]+)[\s]+([^\s]+[^\r\n]*)/){
				$desc = $1;
			}
			${$ret[$#ret]}{"name"} = $name;
			${$ret[$#ret]}{"desc"} = $desc;
			${$ret[$#ret]}{"seq"} = "";
		}else{
			$ss =~ s/[\s]//g;
			${$ret[$#ret]}{"seq"} .= $ss;
		}
	}
	close(IN);
	
	if(length(${$ret[0]}{"seq"}) == 0){
		shift(@ret);
	}
	
	return \@ret;
}

	
