use strict;
use warnings;


require "smithwaterman.pl";
#make formatted output of smithwaterman.pl
my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

my %s1 = %{${loadFasta($file1)}[0]};
my %s2 = %{${loadFasta($file2)}[0]};

my ($res1,$res2,$score) = sw::align($s1{"seq"},$s2{"seq"});
print "#score: $score\n";
print $s1{"name"}." ".join("",@{$res1})."\n";
print $s2{"name"}." ".join("",@{$res2})."\n";

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

	
