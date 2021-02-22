use strict;
use warnings;



#foreach my $ee(keys %answers){
#	print $ee."\n";
#	print ${$answers{$ee}}[0]."\n";
#	print ${$answers{$ee}}[1]."\n";
#}
#exit(0);

my %sourcefas = getSourceFasta();
my %answers = getEmbossResult();

my $okcount = 0;

my $check_pl = 1;
my $check_rs = 1;
my $check_rs_ocl = 1;
my $check_py = 1;
my $check_java = 1;


if($check_pl == 1){
	my $plres  ="res_pl/";
	opendir(DIR,$plres);
	my @allfiles = grep(/.res$/,readdir(DIR));
	closedir(DIR);

	foreach my $aa(@allfiles){
		my $infile = $plres."/".$aa;
		my @r = getSeqOneline($infile);
		my $seq1 = $r[0];
		my $seq2 = $r[1];
		
		
		my $tag = $aa;
		$tag =~ s/\.res.*//;
		my $p1 = $seq1;
		$p1 =~ s/[^A-Z]//g;
		my $p2 = $seq2;
		$p2 =~ s/[^A-Z]//g;
		if($p1 ne ${$sourcefas{$tag}}[0]){
			print $infile."\n";
			print $p1."\n";
			print ${$sourcefas{$tag}}[0]."\n";
			if($p1 !~ /[JUZBOX]/){
				die;
			}
		}
		if($p2 ne ${$sourcefas{$tag}}[1]){
			print $infile."\n";
			print $p2."\n";
			print ${$sourcefas{$tag}}[1]."\n";
			if($p2 !~ /[JUZBOX]/){
				die;
			}
		}
		
		
		my @res = trimTerminal($seq1,$seq2);
		if($aa =~ /seq([0-9]+)/){
			my $fname = "res".$1.".dat";
			
			if($res[0] ne ${$answers{$fname}}[0]){
				print $infile."\n";
				print $res[0]."\n";
				print ${$answers{$fname}}[0]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
					die;
				}
			}
			if($res[1] ne ${$answers{$fname}}[1]){
				print $infile."\n";
				print $res[1]."\n";
				print ${$answers{$fname}}[1]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
					die;
				}
			}
			$okcount ++;
		}else{
			die;
		}
		
	}
}


if($check_py == 1){

	my $pyres  ="res_py/";
	opendir(DIR,$pyres);
	my @allfiles = grep(/.res$/,readdir(DIR));
	closedir(DIR);

	foreach my $aa(@allfiles){
		my $infile = $pyres."/".$aa;
		my @r = getSeqOneline($infile);
		my $seq1 = $r[0];
		my $seq2 = $r[1];
		
		
		my $tag = $aa;
		$tag =~ s/\.res.*//;
		my $p1 = $seq1;
		$p1 =~ s/[^A-Z]//g;
		my $p2 = $seq2;
		$p2 =~ s/[^A-Z]//g;
		if($p1 ne ${$sourcefas{$tag}}[0]){
			print $infile."\n";
			print $p1."\n";
			print ${$sourcefas{$tag}}[0]."\n";
			if($p1 !~ /[JUZBOX]/){
				die;
			}
		}
		if($p2 ne ${$sourcefas{$tag}}[1]){
			print $infile."\n";
			print $p2."\n";
			print ${$sourcefas{$tag}}[1]."\n";
			if($p2 !~ /[JUZBOX]/){
				die;
			}
		}
		
		
		
		
		my @res = trimTerminal($seq1,$seq2);
		
		
		if($aa =~ /seq([0-9]+)/){
			my $fname = "res".$1.".dat";
			
			if($res[0] ne ${$answers{$fname}}[0]){
				print $infile."\n";
				print $res[0]."\n";
				print ${$answers{$fname}}[0]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
					die;
				}
			}
			if($res[1] ne ${$answers{$fname}}[1]){
				print $infile."\n";
				print $res[1]."\n";
				print ${$answers{$fname}}[1]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
					die;
				}
			}
			$okcount ++;
		}else{
			die;
		}
		
		
	}
}
if($check_java == 1){

	my $javares  ="res_java/";
	opendir(DIR,$javares);
	my @allfiles = grep(/.res$/,readdir(DIR));
	closedir(DIR);

	foreach my $aa(@allfiles){
		my $infile = $javares."/".$aa;
		my @r = getSeqFastaLike($infile);
		my $seq1 = $r[0];
		my $seq2 = $r[1];
		
		my $tag = $aa;
		$tag =~ s/\.res.*//;
		my $p1 = $seq1;
		$p1 =~ s/[^A-Z]//g;
		my $p2 = $seq2;
		$p2 =~ s/[^A-Z]//g;
		if($p1 ne ${$sourcefas{$tag}}[0]){
			print $infile."\n";
			print $p1."\n";
			print ${$sourcefas{$tag}}[0]."\n";
			if($p1 !~ /[JUZBOX]/){
				die;
			}
		}
		if($p2 ne ${$sourcefas{$tag}}[1]){
			print $infile."\n";
			print $p2."\n";
			print ${$sourcefas{$tag}}[1]."\n";
			if($p2 !~ /[JUZBOX]/){
				die;
			}
		}
		
		
		my @res = trimTerminal($seq1,$seq2);
		if($aa =~ /seq([0-9]+)/){
			my $fname = "res".$1.".dat";
			
			if($res[0] ne ${$answers{$fname}}[0]){
				print $infile."\n";
				print $res[0]."\n";
				print ${$answers{$fname}}[0]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
					die;
				}
			}
			if($res[1] ne ${$answers{$fname}}[1]){
				print $infile."\n";
				print $res[1]."\n";
				print ${$answers{$fname}}[1]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
					die;
				}
			}
			$okcount ++;
		}else{
			die;
		}
		
	}
}



if($check_rs == 1){

	my $rustres  ="res_rust/";
	opendir(DIR,$rustres);
	my @allfiles = grep(/.res..+$/,readdir(DIR));
	closedir(DIR);

	foreach my $aa(@allfiles){

		my $infile = $rustres."/".$aa;
		my @r = getSeqFastaLike($infile);
		my $seq1 = $r[0];
		my $seq2 = $r[1];
		
		my $tag = $aa;
		$tag =~ s/\.res.*//;
		my $p1 = $seq1;
		$p1 =~ s/[^A-Z]//g;
		my $p2 = $seq2;
		$p2 =~ s/[^A-Z]//g;
		if($p1 ne ${$sourcefas{$tag}}[0]){
			print $infile."\n";
			print $p1."\n";
			print ${$sourcefas{$tag}}[0]."\n";
			if($p1 !~ /[JUZBOX]/){
				die;
			}
		}
		if($p2 ne ${$sourcefas{$tag}}[1]){
			print $infile."\n";
			print $p2."\n";
			print ${$sourcefas{$tag}}[1]."\n";
			if($p2 !~ /[JUZBOX]/){
				die;
			}
		}
		
		if($aa =~ /glocal/){
			#next;
		}
		my @res;
		if($aa =~ /global/ || $aa =~ /glocal/){
			@res = ($seq1,$seq2);
		}else{
			@res = trimTerminal($seq1,$seq2);
		}
		if($aa =~ /seq([0-9]+)/){
			my $sid = $1;
			my $fname = "res".$sid.".dat";
			if($aa =~ /\.glocal/){
				 $fname = "needle_glocal_res".$sid.".dat";
			}
			if($aa =~ /\.global/){
				 $fname = "needle_res".$sid.".dat";
			}
			
			
			if($res[0] ne ${$answers{$fname}}[0]){
				print $infile."\n";
				print $res[0]."\n";
				print ${$answers{$fname}}[0]."\n";
				print ${$answers{$fname}}[1]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
					die;
				}
			}
			if($res[1] ne ${$answers{$fname}}[1]){
				print $infile."\n";
				print $res[1]."\n";
				print ${$answers{$fname}}[0]."\n";
				print ${$answers{$fname}}[1]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
					die;
				}
			}
			$okcount ++;
		}else{
			die;
		}
		
	}

}


if($check_rs_ocl == 1){

	my $rustres  ="res_rust_ocl/";
	opendir(DIR,$rustres);
	my @allfiles = grep(/.res..+$/,readdir(DIR));
	closedir(DIR);

	foreach my $aa(@allfiles){

		my $infile = $rustres."/".$aa;
		my @r = getSeqFastaLike($infile);
		my $seq1 = $r[0];
		my $seq2 = $r[1];
		
		my $tag = $aa;
		$tag =~ s/\.res.*//;
		my $p1 = $seq1;
		$p1 =~ s/[^A-Z]//g;
		my $p2 = $seq2;
		$p2 =~ s/[^A-Z]//g;
		if($p1 ne ${$sourcefas{$tag}}[0]){
			print $infile."\n";
			print $p1."\n";
			print ${$sourcefas{$tag}}[0]."\n";
			if($p1 !~ /[JUZBOX]/){
				die;
			}
		}
		if($p2 ne ${$sourcefas{$tag}}[1]){
			print $infile."\n";
			print $p2."\n";
			print ${$sourcefas{$tag}}[1]."\n";
			if($p2 !~ /[JUZBOX]/){
				die;
			}
		}
		
		if($aa =~ /glocal/){
			#next;
		}
		my @res;
		if($aa =~ /global/ || $aa =~ /glocal/){
			@res = ($seq1,$seq2);
		}else{
			@res = trimTerminal($seq1,$seq2);
		}
		if($aa =~ /seq([0-9]+)/){
			my $sid = $1;
			my $fname = "res".$sid.".dat";
			if($aa =~ /\.glocal/){
				 $fname = "needle_glocal_res".$sid.".dat";
			}
			if($aa =~ /\.global/){
				 $fname = "needle_res".$sid.".dat";
			}
			
			
			if($res[0] ne ${$answers{$fname}}[0]){
				print $infile."\n";
				print $res[0]."\n";
				print ${$answers{$fname}}[0]."\n";
				print ${$answers{$fname}}[1]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
				#	die;
				}
			}
			if($res[1] ne ${$answers{$fname}}[1]){
				print $infile."\n";
				print $res[1]."\n";
				print ${$answers{$fname}}[0]."\n";
				print ${$answers{$fname}}[1]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
				#	die;
				}
			}
			$okcount ++;
		}else{
			die;
		}
		
	}
}

if($check_rs_ocl == 1){

	my $rustres  ="res_rust_ocl_list/";
	opendir(DIR,$rustres);
	my @allfiles = grep(/.res..+$/,readdir(DIR));
	closedir(DIR);

	foreach my $aa(@allfiles){

		my $infile = $rustres."/".$aa;
		my @r = getSeqFastaLike($infile);
		my $seq1 = $r[0];
		my $seq2 = $r[1];
		
		my $tag = $aa;
		$tag =~ s/\.res.*//;
		my $p1 = $seq1;
		$p1 =~ s/[^A-Z]//g;
		my $p2 = $seq2;
		$p2 =~ s/[^A-Z]//g;
		if($p1 ne ${$sourcefas{$tag}}[0]){
			print $infile."\n";
			print $p1."\n";
			print ${$sourcefas{$tag}}[0]."\n";
			if($p1 !~ /[JUZBOX]/){
				die;
			}
		}
		if($p2 ne ${$sourcefas{$tag}}[1]){
			print $infile."\n";
			print $p2."\n";
			print ${$sourcefas{$tag}}[1]."\n";
			if($p2 !~ /[JUZBOX]/){
				die;
			}
		}
		
		if($aa =~ /glocal/){
			#next;
		}
		my @res;
		if($aa =~ /global/ || $aa =~ /glocal/){
			@res = ($seq1,$seq2);
		}else{
			@res = trimTerminal($seq1,$seq2);
		}
		if($aa =~ /seq([0-9]+)/){
			my $sid = $1;
			my $fname = "res".$sid.".dat";
			if($aa =~ /\.glocal/){
				 $fname = "needle_glocal_res".$sid.".dat";
			}
			if($aa =~ /\.global/){
				 $fname = "needle_res".$sid.".dat";
			}
			
			
			if($res[0] ne ${$answers{$fname}}[0]){
				print $infile."\n";
				print $res[0]."\n";
				print ${$answers{$fname}}[0]."\n";
				print ${$answers{$fname}}[1]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
				#	die;
				}
			}
			if($res[1] ne ${$answers{$fname}}[1]){
				print $infile."\n";
				print $res[1]."\n";
				print ${$answers{$fname}}[0]."\n";
				print ${$answers{$fname}}[1]."\n";
				if($res[0] !~ /[JUZBOX]/ && $res[1] !~ /[JUZBOX]/){
				#	die;
				}
			}
			$okcount ++;
		}else{
			die;
		}
		
	}
}

print "OK\nChecked $okcount results.\n";






sub trimTerminal{
	my @s1 = split(//,$_[0]);
	my @s2 = split(//,$_[1]);
	for(my $ii = 0;$ii <= $#s1;$ii++){
		if($s1[$ii] ne "-" && $s2[$ii] ne "-"){
			last;
		}
		$s1[$ii] = "";
		$s2[$ii] = "";
	}
	
	for(my $i = 0;$i <= $#s1;$i++){
		my $ii = $#s1 -$i;
		if($s1[$ii] ne "-" && $s2[$ii] ne "-"){
			last;
		}
		$s1[$ii] = "";
		$s2[$ii] = "";
	}
	
	
	return (join("",@s1),join("",@s2));
}


sub getEmbossResult{
	my $embossres  ="emboss_results/";
	opendir(DIR,$embossres);
	my @emm = grep(/.dat$/,readdir(DIR));
	closedir(DIR);
	my %answers;
	foreach my $ee(@emm){
		my @tmpp;
		my $seq1 = "";
		my $seq2 = "";
		open(IN,$embossres."/".$ee);
		while(my $ss = <IN>){
			if($ss =~ /^[\s]*s1[\s]+[0-9]+[\s]*([^\s]+)/){
				$seq1 .= $1;
			}
			if($ss =~ /^[\s]*s2[\s]+[0-9]+[\s]*([^\s]+)/){
				$seq2 .= $1;
			}
		}
		close(IN);
		push(@tmpp,$seq1);
		push(@tmpp,$seq2);
		$answers{$ee} = \@tmpp;
	}
	return %answers;
}

sub getSourceFasta{
	my %ret;
	my $seqdir  ="inputs/";
	opendir(DIR,$seqdir);
	my @allfiles = grep(/\.1\.fas$/,readdir(DIR));
	closedir(DIR);
	foreach my $ss(@allfiles){
		my $infile1 = $seqdir."/".$ss;
		my $infile2 = $seqdir."/".$ss;
		$infile2 =~ s/\.1\.fas/.2.fas/;
		my $tag = $ss;
		my @s1 = getSeqFastaLike($infile1);
		my @s2 = getSeqFastaLike($infile2);
		my @tmp;
		push(@tmp,$s1[0]);
		push(@tmp,$s2[1]);
		$ret{$tag} = \@tmp;
	}
	return %ret;
}



sub getSeqFastaLike{
	my $filename = $_[0];
	my $flag = 0;
	my $s1 = "";
	my $s2 = "";
	open(IN,$filename);
	while(my $ss = <IN>){
		if($ss =~ />s1/){
			$flag = 1;
			next;
		}
		if($ss =~ />s2/){
			$flag = 2;
			next;
		}
		if($flag == 1){
			$s1 .= $ss;
		}
		if($flag == 2){
			$s2 .= $ss;
		}
	}
	close(IN);
	
	$s1 =~ s/[\s]//g;
	$s2 =~ s/[\s]//g;
	return ($s1,$s2);
}

sub getSeqOneline{
	my $filename = $_[0];
	open(IN,$filename);
	my $s1 = "";
	my $s2 = "";
	while(my $ss = <IN>){
		if($ss =~ /^[\s]*s1[\s]+([^\s]+)/){
			$s1 .= $1;
		}
		if($ss =~ /^[\s]*s2[\s]+([^\s]+)/){
			$s2 .= $1;
		}
	}
	close(IN);
	$s1 =~ s/[\s]//g;
	$s2 =~ s/[\s]//g;
	return ($s1,$s2);
}