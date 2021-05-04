extern crate sa_opencl;

use std::{env};
use sa_opencl::opencl_sequence_alignment;
use std::io::{BufReader,BufRead,BufWriter,Write};
use std::fs::File;


//https://stackoverflow.com/questions/37888042/remove-single-trailing-newline-from-string-without-cloning
//Sven Marnach

fn trim_newline(s: &mut String) {
    if s.ends_with('\n') {
        s.pop();
        if s.ends_with('\r') {
            s.pop();
        }
    }
}

struct AlignmentOptions{
    alignment_type:usize,
    file1:String,
    file2:String,
    outfilename:String,
    list:bool,
    clustering:bool,
    c_identity:Option<f64>,
    c_coverage_short:Option<f64>,
    c_coverage_long:Option<f64>,

}

impl AlignmentOptions{
    fn parse(argss:&Vec<String>)->AlignmentOptions{
        let mut ret = AlignmentOptions{
            alignment_type:opencl_sequence_alignment::ALIGN_LOCAL,
            file1:"".to_owned(),
            file2:"".to_owned(),
            list:false,
            clustering:false,
            c_identity:None,
            c_coverage_short:None,
            c_coverage_long:None,
            outfilename:"".to_owned()
        };
        let mut flag:Vec<bool> = vec![false;argss.len()];
        let mut file_candidates:Vec<String> = vec![];
        for ii in 0..argss.len(){
            if argss[ii] == "-glocal" || argss[ii] == "-global" || argss[ii] == "-local"{
                if argss[ii] == "-glocal"{
                    ret.alignment_type = opencl_sequence_alignment::ALIGN_GLOCAL;
                }else if argss[ii] == "-global"{
                    ret.alignment_type = opencl_sequence_alignment::ALIGN_GLOBAL;
                }else{
                    ret.alignment_type = opencl_sequence_alignment::ALIGN_LOCAL;
                }
                flag[ii] = true;
            }else if argss[ii] == "-list"{
                ret.list = true;
                flag[ii] = true;
            }else if argss[ii] == "-cluster" || argss[ii] == "-clustering" {
                ret.clustering = true;
                flag[ii] = true;
            }else if argss[ii] == "-coverage_short"{
                ret.c_coverage_short = Some(argss[ii+1].parse::<f64>().unwrap_or_else(|e|panic!("parse error {} {:?}",argss[ii+1],e)));
                flag[ii] = true;
                flag[ii+1] = true;
            }else if argss[ii] == "-coverage_long"{
                ret.c_coverage_long =  Some(argss[ii+1].parse::<f64>().unwrap_or_else(|e|panic!("parse error {} {:?}",argss[ii+1],e)));
                flag[ii] = true;
                flag[ii+1] = true;
            }else if argss[ii] == "-identity"{
                flag[ii] = true;
                flag[ii+1] = true;
                ret.c_identity = Some(argss[ii+1].parse::<f64>().unwrap_or_else(|e|panic!("parse error {} {:?}",argss[ii+1],e)));
            }else if argss[ii] == "-out"{
                flag[ii] = true;
                flag[ii+1] = true;
                ret.outfilename = argss[ii+1].to_string();
            }else{
                if !flag[ii]{
                    if argss[ii].get(0..1).unwrap() == "-"{
                        panic!("Unknown option {}",argss[ii]);
                    }else{
                        file_candidates.push(argss[ii].to_string());
                    }
                }
            }
        }
        if !ret.clustering && !ret.list{
            if file_candidates.len() != 2{
                panic!("2 files must be provided {:?}.",file_candidates);
            }
            ret.file1 = file_candidates[0].clone();
            ret.file2 = file_candidates[1].clone();
        }else{
            if ret.clustering && ret.list{
                panic!("Incompatible option -list & -cluster(ing)");
            }
            if file_candidates.len() != 1{
                panic!("1 file must be provided {:?}.",file_candidates);
            }
            ret.file1 = file_candidates[0].clone();
            if ret.clustering{
                if ret.outfilename.len() == 0{
                    panic!("Clustering must have -out.");
                }
            }
        }
        return ret;
    }

}

fn main(){
	let argss: Vec<String> = env::args().collect();
    let mut mess:&str = "Local alignment";
    if argss.len() < 2{
        eprintln!("usage: sa_opencl [(-global|-glocal|-local(default))] <infile1 (fasta file)>  <infile2 (fasta file)> ");
        eprintln!("usage: sa_opencl [(-global|-glocal|-local(default))] [-list] <list file>");
        eprintln!("usage: sa_opencl -cluster [-identity 0.0-1.0] [-coverage_short 0.0-1.0] [-coverage_long 0.0-1.0] <fasta file>");
        eprintln!("The \"list file\" has a list of tab separated pairs as follows.");
        eprintln!("<infile1 (fasta file)>  <infile2 (fasta file)>");
        eprintln!("<infile3 (fasta file)>  <infile4 (fasta file)>");
        eprintln!("<infile5 (fasta file)>  <infile6 (fasta file)>");
        eprintln!("...");
        eprintln!("Then,");
        eprintln!("sequences in infile1 and infile2,");
        eprintln!("sequences in infile3 and infile4,");
        eprintln!("sequences in infile5 and infile6,");
        eprintln!("will be aligned.");
        std::process::exit(-1);
    }
    
    let opp = AlignmentOptions::parse(&argss);
    if opp.alignment_type == opencl_sequence_alignment::ALIGN_GLOCAL{
        mess = "Glocal alignment";
    }else if opp.alignment_type == opencl_sequence_alignment::ALIGN_GLOBAL{
        mess = "Global alignment";
    }else if opp.alignment_type == opencl_sequence_alignment::ALIGN_LOCAL{
        mess = "Local alignment";
    }
    let mut sw = opencl_sequence_alignment::OpenCLSequenceAlignment::new(1000,Box::new(opencl_sequence_alignment::SubstitutionMatrix::get_blosum62_matrix()),10.0,0.5,opp.alignment_type);
    
    let mut filelist:Vec<(String,String)> = vec![];
    if opp.clustering{
        let mut seq1 = opencl_sequence_alignment::SeqData::load_fasta(&opp.file1,false);
        seq1.sort_by(|a,b|b.seq.len().cmp(&a.seq.len()));
        let mut cluster_of:Vec<usize> = (0..seq1.len()).into_iter().collect();
        let mut members:Vec<Vec<usize>> = vec![vec![];cluster_of.len()];
        let coverage_short = opp.c_coverage_short.unwrap_or(0.8);
        let coverage_long = opp.c_coverage_long.unwrap_or(0.8);
        let identity_threshold:f64 = opp.c_identity.unwrap_or(0.8);

        for ii in 0..seq1.len(){
            if cluster_of[ii] != ii{
                continue;
            }
            members[ii].push(ii);
            for jj in (ii+1)..seq1.len(){
                if cluster_of[jj] != jj{
                    continue;
                }
                let res = sw.align(&seq1[ii],&seq1[jj],false);
                let mut alen:usize = 0;
                let mut blen:usize = 0;
                let mut matchnum:usize = 0;
                for kk in 0..res.0.len(){
                    if res.0[kk] != "-"{
                        alen += 1;
                    }
                    if res.1[kk] != "-"{
                        blen += 1;
                    }
                    if res.0[kk] == "-" || res.1[kk] == "-"{
                        continue;
                    }
                    if res.0[kk] == res.1[kk]{
                        matchnum += 1;
                    }
                }
                let lcov:f64;
                let scov:f64;
                if seq1[ii].seq.len() > seq1[jj].seq.len(){
                    lcov = alen as f64/(seq1[ii].seq.len() as f64);
                    scov = blen as f64/(seq1[jj].seq.len() as f64);
                }else{
                    scov = alen as f64/(seq1[ii].seq.len() as f64);
                    lcov = blen as f64/(seq1[jj].seq.len() as f64);
                }
                let ident = matchnum as f64/(res.0.len() as f64);
                if lcov >= coverage_long && scov >= coverage_short && ident >= identity_threshold{
                    cluster_of[jj] = ii;
                    members[ii].push(jj);
                }
            }
        }
        let mut f = BufWriter::new(File::create(opp.outfilename.clone()).unwrap());
        for cc in 0..cluster_of.len(){
            if cc == cluster_of[cc]{
                let s1 = seq1[cc].seq.iter().fold("".to_string(),|s,m|s+m);
                f.write(format!(">{} {}\n{}\n",seq1[cc].name,seq1[cc].desc,s1).as_bytes()).unwrap(); 
            }
        }
        f.flush().unwrap();
        let mut ff = BufWriter::new(File::create(opp.outfilename.clone()+".clstr").unwrap());
        for cc in 0..cluster_of.len(){
            for (mii,mm) in members[cc].iter().enumerate(){
                if mii != 0{
                    ff.write(" ".as_bytes()).unwrap();
                }
                ff.write(format!("{}",seq1[*mm].name).as_bytes()).unwrap();
            }
            ff.write("\n".as_bytes()).unwrap();
        }
        ff.flush().unwrap();
        

    }else{
        if opp.list{
            let file = File::open(opp.file1).unwrap();
            let reader = BufReader::new(file);
            for (_lcount,line_) in reader.lines().enumerate() {
                let mut line =line_.unwrap();
                trim_newline(&mut line);
                let mut spp:Vec<String> = line.as_str().split("\t").into_iter().map(|m|m.to_owned()).collect();
                if spp.len() == 1{
                    spp = line.as_str().split(" ").into_iter().map(|m|m.to_owned()).collect();
                }
                if spp.len() > 2{
                    println!("{} \n^ Only {} {} are used.",line,spp[0],spp[1]);
                }else{
                    if spp.len() < 2{
                        println!("{} \n is ignoed.",line);
                        continue;
                    }
                    filelist.push((spp[0].clone(),spp[1].clone()));
                }
            }
        }else{
            filelist.push((opp.file1.clone(),opp.file2.clone()));    
        }
        let mut outfile:Option<BufWriter<File>> = if opp.outfilename.len() > 0{
            Some(BufWriter::new(File::create(opp.outfilename).unwrap()))
        }else{
            None
        };
        for (file1,file2) in filelist.into_iter(){
            let seq1 = opencl_sequence_alignment::SeqData::load_fasta(&file1,false);
            let seq2 = opencl_sequence_alignment::SeqData::load_fasta(&file2,false);
            for ss1 in seq1.iter(){
                for ss2 in seq2.iter(){
                    let res = sw.align(ss1,ss2,true);
                    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
                    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
                    if let Some(f) = outfile.as_mut(){
                        f.write_all(format!("#score:{}",res.2).as_bytes()).unwrap();
                        f.write_all(format!("#type:{}",mess).as_bytes()).unwrap();
                        f.write_all(format!(">{}\n{}\n",ss1.name,r1).as_bytes()).unwrap();
                        f.write_all(format!(">{}\n{}\n",ss2.name, r2).as_bytes()).unwrap();
                    }else{
                        println!("#score:{}",res.2);
                        println!("#type:{}",mess);
                        println!(">{}\n{}\n",ss1.name,r1);
                        println!(">{}\n{}\n",ss2.name, r2);
                    }
                }
            }
        }
    }
    return;
}