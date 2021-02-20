extern crate sa_opencl;

use std::env;
use sa_opencl::opencl_sequence_alignment;
use std::io::{BufReader,BufRead};
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

fn main(){
	let argss: Vec<String> = env::args().collect();
    let mut alignment_type:usize = opencl_sequence_alignment::ALIGN_LOCAL;
    let mut mess:&str = "Local alignment";
    if argss.len() < 2{
        eprintln!("usage: sa_opencl [(-global|-glocal|-local(default))] <infile1 (fasta file)>  <infile2 (fasta file)> ");
        eprintln!("usage: sa_opencl [(-global|-glocal|-local(default))] [-list] <list file>");
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
    let mut file1_:&str = &argss[1];
    let mut file2_:&str = &argss[2];

    if file2_ == "-local" || file2_ == "-glocal" || file2_ == "-global"{
        if file1_ == "-list"{
            file1_ = file2_;
            file2_ = "-list";
        }
    }else if file1_ == "-list"{
        file1_ = "-local";
        file2_ = "-list";
    }

    if file1_ == "-global"{
        alignment_type = opencl_sequence_alignment::ALIGN_GLOBAL;
        mess = "Global alignment";
        file1_ = file2_;
        if argss.len() > 3{
            file2_ = &argss[3];
        }else{
            file2_ = &argss[2];
        }
    }else if file1_ == "-glocal"{
        mess = "Glocal alignment";
        alignment_type = opencl_sequence_alignment::ALIGN_GLOCAL;
        file1_ = file2_;
        if argss.len() > 3{
            file2_ = &argss[3];
        }else{
            file2_ = &argss[2];
        }
    }else if file1_ == "-local"{
        mess = "Local alignment";
        alignment_type = opencl_sequence_alignment::ALIGN_LOCAL;
        file1_ = file2_;
        if argss.len() > 3{
            file2_ = &argss[3];
        }else{
            file2_ = &argss[2];
        }
    }

    let mut sw = opencl_sequence_alignment::OpenCLSequenceAlignment::new(1000,Box::new(opencl_sequence_alignment::SubstitutionMatrix::get_blosum62_matrix()),10.0,0.5,alignment_type);
        
    let mut filelist:Vec<(String,String)> = vec![];
    if file1_ == "-list"{

        let file = File::open(file2_).unwrap();
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
        filelist.push((file1_.to_owned(),file2_.to_owned()));    
    }
    for (file1,file2) in filelist.into_iter(){

        let seq1 = opencl_sequence_alignment::SeqData::load_fasta(&file1,false);
        let seq2 = opencl_sequence_alignment::SeqData::load_fasta(&file2,false);
        for ss1 in seq1.iter(){
            for ss2 in seq2.iter(){
                let res = sw.align(ss1,ss2,true);
                let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
                let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);cd 
                println!("#score:{}",res.2);
                println!("#type:{}",mess);
                println!(">{}\n{}\n",ss1.name,r1);
                println!(">{}\n{}\n",ss2.name, r2);            
            }
        }
    }
    return;
}