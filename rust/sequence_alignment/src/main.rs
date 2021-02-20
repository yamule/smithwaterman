extern crate seqalign;

use std::env;
use seqalign::sequence_alignment;

fn main(){
	let argss: Vec<String> = env::args().collect();
    let mut alignment_type:usize = sequence_alignment::ALIGN_LOCAL;
    let mut mess:&str = "Local alignment";
    if argss.len() < 2{
        panic!("usage: sequence_alignment <infile1 (single plain fasta file)>  <infile2 (single plain fasta file)> ");
    }
    let mut file1:&str = &argss[1];
    let mut file2:&str = &argss[2];
    
    if file1 == "-global"{
        alignment_type = sequence_alignment::ALIGN_GLOBAL;
        mess = "Global alignment";
        file1 = &argss[2];
        file2 = &argss[3];
    }else if file1 == "-glocal"{
        mess = "Glocal alignment";
        alignment_type = sequence_alignment::ALIGN_GLOCAL;
        file1 = &argss[2];
        file2 = &argss[3];
    }else if file1 == "-local"{
        mess = "Local alignment";
        alignment_type = sequence_alignment::ALIGN_LOCAL;
        file1 = &argss[2];
        file2 = &argss[3];
    }
    let seq1 = sequence_alignment::SeqData::load_fasta(file1,false);
    let seq2 = sequence_alignment::SeqData::load_fasta(file2,false);
    let mut sw = sequence_alignment::SequenceAlignment::new(Box::new(sequence_alignment::SubstitutionMatrix::get_blosum62_matrix()),10.0,0.5,alignment_type);
    for ss1 in seq1.iter(){
        for ss2 in seq2.iter(){
            let res = sw.align(ss1,ss2,true);
            let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
            let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
            println!("#score:{}",res.2);
            println!("#type:{}",mess);
            println!(">{}\n{}\n",ss1.name,r1);
            println!(">{}\n{}\n",ss2.name, r2);            
        }
    }
    return;
}
