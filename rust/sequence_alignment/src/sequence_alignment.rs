
use std::fs::File;
use std::collections::HashMap;

use std::io::{BufReader,BufRead};

pub const CELL_MATCH:usize = 0;
pub const CELL_GAPINX:usize = 1;
pub const CELL_GAPINY:usize = 2;

pub const ALIGN_GLOBAL:usize = 0;
pub const ALIGN_GLOCAL:usize = 1;
pub const ALIGN_LOCAL:usize = 2;

pub struct SequenceAlignment{
    pub cells:Vec<Vec<SWCell>>,
    pub scoring_matrix:Box<dyn ScoringMatrix>,
    pub e_penalty:f32,
    pub o_penalty:f32,
    pub seqlen_a:usize,
    pub seqlen_b:usize,
    pub alignment_type:usize,
}

impl SequenceAlignment{
    pub fn new(ss:Box<dyn ScoringMatrix>,go:f32,ge:f32,alignment_type:usize)->SequenceAlignment{

        let mut pgo = go;
        let mut pge = ge;
        if pgo > 0.0{
            pgo *= -1.0_f32;
        }
        if pge > 0.0{
            pge *= -1.0_f32;
        }
        //abs()*-1.0 とかにしてもいいと思うが、丸め誤差が出ると嫌なので・・・

        return SequenceAlignment{cells:vec![vec![]],scoring_matrix:ss
            ,o_penalty:pgo
            ,e_penalty:pge
            ,seqlen_a:0
            ,seqlen_b:0
            ,alignment_type:alignment_type
        };
    }

    pub fn clear_buffer(&mut self){
        self.cells = vec![vec![]];
    }

    pub fn prepare(&mut self,s1:&SeqData,s2:&SeqData){
        //substitution matrix の場合は必要ない
        self.scoring_matrix.prepare(&s1,&s2);
    }
    pub fn fill_matrix(&mut self,seq1:&Vec<usize>,seq2:&Vec<usize>
    ,partial_region:Option<(usize,usize)>){
        if let Some(x) = partial_region{
            self.seqlen_a = x.0;
            self.seqlen_b = x.1;
        }else{
            let len1 = seq1.len();
            let len2 = seq2.len();
            self.seqlen_a = len1;
            self.seqlen_b = len2;
        }
        
        if self.cells.len() == 0{
            self.cells = vec![vec![SWCell::new();self.seqlen_b+1];self.seqlen_a+1];
        }else if self.seqlen_a+1 > self.cells.len() || self.seqlen_b+1 > self.cells[0].len(){
            self.cells = vec![vec![SWCell::new();self.seqlen_b+1];self.seqlen_a+1];
        }

        //テストしてない
        let start_epenal:f32;
        let start_openal:f32;
        if self.alignment_type == ALIGN_GLOBAL{
            start_epenal = self.e_penalty;
            start_openal = self.o_penalty;
        }else{
            start_epenal = 0.0_f32;
            start_openal = 0.0_f32;
        }

        for ii in 0..(self.seqlen_a+1){    
            for jj in 0..(self.seqlen_b+1){
                let mut openal = self.o_penalty;
                let mut epenal = self.e_penalty;
                if ii == 0 && jj == 0{
                    self.cells[ii][jj].set(CELL_MATCH,CELL_MATCH,0.0);
                    self.cells[ii][jj].set(CELL_GAPINX,CELL_MATCH,-1.0);
                    self.cells[ii][jj].set(CELL_GAPINY,CELL_MATCH,-1.0);
                    continue;
                }

                //if ii*jj*(self.seqlen_a-ii)*(self.seqlen_b-jj) == 0{
                if ii*jj == 0{
                    openal = start_openal;
                    epenal = start_epenal;
                }
                if ii == 0{
                    let lscore = jj as f32*epenal+(openal-epenal);
                    self.cells[ii][jj].scores[CELL_MATCH ] = lscore+10.0*self.o_penalty +10.0*self.e_penalty ;
                    self.cells[ii][jj].scores[CELL_GAPINX] = lscore;
                    self.cells[ii][jj].scores[CELL_GAPINY] = lscore+10.0*self.o_penalty +10.0*self.e_penalty ;
                    self.cells[ii][jj].prev[CELL_MATCH ] = CELL_GAPINX;
                    self.cells[ii][jj].prev[CELL_GAPINX] = CELL_GAPINX;
                    self.cells[ii][jj].prev[CELL_GAPINY] = CELL_GAPINX;
                    continue;
                }else if jj == 0{
                    let lscore = ii as f32*epenal+(openal-epenal);
                    self.cells[ii][jj].scores[CELL_MATCH ] = lscore+10.0*self.o_penalty +10.0*self.e_penalty ;
                    self.cells[ii][jj].scores[CELL_GAPINX] = lscore+10.0*self.o_penalty +10.0*self.e_penalty ;
                    self.cells[ii][jj].scores[CELL_GAPINY] = lscore;
                    self.cells[ii][jj].prev[CELL_MATCH ] = CELL_GAPINY;
                    self.cells[ii][jj].prev[CELL_GAPINX] = CELL_GAPINY;
                    self.cells[ii][jj].prev[CELL_GAPINY] = CELL_GAPINY;
                    continue;
                }


                let cel_lt = &self.cells[ii-1][jj-1];
                let cel_l = &self.cells[ii-1][jj];
                let cel_t = &self.cells[ii][jj-1];
                


                let mmscore:f32 = self.scoring_matrix.get_score(seq1[ii-1],seq2[jj-1]);
                
                let mut matchscore:f32;
                let matchindex:usize;

                let mut gapxscore:f32;
                let gapxindex:usize;

                let mut gapyscore:f32;
                let gapyindex:usize;

                //Cell にマッチした場合
                if  cel_lt.scores[CELL_MATCH]+mmscore
                >= cel_lt.scores[CELL_GAPINX]+mmscore{
                    if cel_lt.scores[CELL_MATCH]+mmscore
                    >= cel_lt.scores[CELL_GAPINY]+mmscore{    
                        matchscore = cel_lt.scores[CELL_MATCH]+mmscore;
                        matchindex = CELL_MATCH;
                    }else{
                        matchscore = cel_lt.scores[CELL_GAPINY]+mmscore;
                        matchindex = CELL_GAPINY;
                    }
                }else{
                    if cel_lt.scores[CELL_GAPINX]+mmscore
                    >= cel_lt.scores[CELL_GAPINY]+mmscore{    
                        matchscore = cel_lt.scores[CELL_GAPINX]+mmscore;
                        matchindex = CELL_GAPINX;
                    }else{
                        matchscore = cel_lt.scores[CELL_GAPINY]+mmscore;
                        matchindex = CELL_GAPINY;
                    }
                }





                if self.alignment_type != ALIGN_LOCAL{
                    //emboss の water と needle で優先順位が違うようだ・・・
                    let popenal = if self.seqlen_a-ii == 0{start_openal}else{openal};
                    let pepenal = if self.seqlen_a-ii == 0{start_epenal}else{epenal};
                    
                    let qopenal = if self.seqlen_b-jj == 0{start_openal}else{openal};
                    let qepenal = if self.seqlen_b-jj == 0{start_epenal}else{epenal};
                    
                    if  cel_t.scores[CELL_MATCH]+popenal
                    > cel_t.scores[CELL_GAPINX]+pepenal{
                        if cel_t.scores[CELL_MATCH]+popenal >= cel_t.scores[CELL_GAPINY]+popenal{
                            gapxscore = cel_t.scores[CELL_MATCH]+popenal;
                            gapxindex = CELL_MATCH;
                        }else{
                            gapxscore = cel_t.scores[CELL_GAPINY]+popenal;
                            gapxindex = CELL_GAPINY;
                        }
                    }else{
                        if cel_t.scores[CELL_GAPINX]+pepenal >= cel_t.scores[CELL_GAPINY]+popenal{
                            gapxscore = cel_t.scores[CELL_GAPINX]+pepenal;
                            gapxindex = CELL_GAPINX;
                        }else{
                            gapxscore = cel_t.scores[CELL_GAPINY]+popenal;
                            gapxindex = CELL_GAPINY;
                        }
                    }


                    if  cel_l.scores[CELL_MATCH]+qopenal
                    > cel_l.scores[CELL_GAPINY]+qepenal{
                        if cel_l.scores[CELL_MATCH]+qopenal >= cel_l.scores[CELL_GAPINX]+qopenal{
                            gapyscore = cel_l.scores[CELL_MATCH]+qopenal;
                            gapyindex = CELL_MATCH;
                        }else{
                            gapyscore = cel_l.scores[CELL_GAPINX]+qopenal;
                            gapyindex = CELL_GAPINX;
                            
                        }

                    }else{
                        if cel_l.scores[CELL_GAPINY]+qepenal >= cel_l.scores[CELL_GAPINX]+qopenal{
                            gapyscore = cel_l.scores[CELL_GAPINY]+qepenal;
                            gapyindex = CELL_GAPINY;
                        }else{
                            gapyscore = cel_l.scores[CELL_GAPINX]+qopenal;
                            gapyindex = CELL_GAPINX;
                        }
                    }
                    
                }else{
                    if  cel_t.scores[CELL_MATCH]+openal
                    >= cel_t.scores[CELL_GAPINX]+epenal{
                        if cel_t.scores[CELL_MATCH]+openal > cel_t.scores[CELL_GAPINY]+openal{
                            gapxscore = cel_t.scores[CELL_MATCH]+openal;
                            gapxindex = CELL_MATCH;
                        }else{
                            gapxscore = cel_t.scores[CELL_GAPINY]+openal;
                            gapxindex = CELL_GAPINY;
                        }
                    }else{
                        if cel_t.scores[CELL_GAPINX]+epenal > cel_t.scores[CELL_GAPINY]+openal{
                            gapxscore = cel_t.scores[CELL_GAPINX]+epenal;
                            gapxindex = CELL_GAPINX;
                        }else{
                            gapxscore = cel_t.scores[CELL_GAPINY]+openal;
                            gapxindex = CELL_GAPINY;
                        }
                    }

                    if  cel_l.scores[CELL_MATCH]+openal
                    >= cel_l.scores[CELL_GAPINY]+epenal{
                        if cel_l.scores[CELL_MATCH]+openal > cel_l.scores[CELL_GAPINX]+openal{
                            gapyscore = cel_l.scores[CELL_MATCH]+openal;
                            gapyindex = CELL_MATCH;
                        }else{
                            gapyscore = cel_l.scores[CELL_GAPINX]+openal;
                            gapyindex = CELL_GAPINX;
                            
                        }

                    }else{
                        if cel_l.scores[CELL_GAPINY]+epenal > cel_l.scores[CELL_GAPINX]+openal{
                            gapyscore = cel_l.scores[CELL_GAPINY]+epenal;
                            gapyindex = CELL_GAPINY;
                        }else{
                            gapyscore = cel_l.scores[CELL_GAPINX]+openal;
                            gapyindex = CELL_GAPINX;
                        }
                    }
                }
                
                //可読性のためであり上のブロックに入れてもよいと思う
                if self.alignment_type == ALIGN_LOCAL{
                    matchscore = matchscore.max(0.0);
                    gapxscore = gapxscore.max(0.0);
                    gapyscore = gapyscore.max(0.0);
                }

                self.cells[ii][jj].scores[CELL_MATCH ] = matchscore;
                self.cells[ii][jj].scores[CELL_GAPINX] = gapxscore;
                self.cells[ii][jj].scores[CELL_GAPINY] = gapyscore;
                
                self.cells[ii][jj].prev[CELL_MATCH ] = matchindex;
                self.cells[ii][jj].prev[CELL_GAPINX] = gapxindex;
                self.cells[ii][jj].prev[CELL_GAPINY] = gapyindex;
            }    
        }

    }
    pub fn backtrack(&mut self)->(Vec<i64>,Vec<i64>,f32){
        let mut startx_:i64 = -1;
        let mut starty_:i64 = -1;
        let mut maxscore:f32;
        let mut max_place:usize;

        let mut ret1:Vec<i64> = vec![];
        let mut ret2:Vec<i64> = vec![];

        if self.alignment_type == ALIGN_LOCAL{
            let xlen = self.seqlen_a +1;
            let ylen = self.seqlen_b +1;
            maxscore = 0.0;
            for ii in 0..xlen{
                for jj in 0..ylen{
                    if self.cells[ii][jj].scores[CELL_MATCH] > maxscore{
                        maxscore = self.cells[ii][jj].scores[CELL_MATCH];
                        startx_ = ii as i64;
                        starty_ = jj as i64;
                    }
                }
            }
            max_place = CELL_MATCH;
        }else if self.alignment_type == ALIGN_GLOBAL || self.alignment_type == ALIGN_GLOCAL{
            let xlen = self.seqlen_a +1;
            let ylen = self.seqlen_b +1;
            startx_ = xlen as i64-1;
            starty_ = ylen as i64-1;
            let slen = self.cells[startx_ as usize][starty_ as usize].scores.len();
            maxscore = self.cells[startx_ as usize][starty_ as usize].scores[0];
            max_place = 0;
            for ii in 1..slen{
                if maxscore < self.cells[startx_ as usize][starty_ as usize].scores[ii]{
                    maxscore = self.cells[startx_ as usize][starty_ as usize].scores[ii];
                    max_place = ii;
                }
            }
        }else{
            //glocal に使おうと思っていたが今使ってない
            panic!("not used");
            /*
            let xlen = self.seqlen_a +1;
            let ylen = self.seqlen_b +1;
            assert_eq!(3,self.cells[xlen-1][0].scores.len());
            maxscore = self.cells[xlen-1][0].scores[2];
            max_place = CELL_MATCH;
            startx_ = xlen as i64-1;
            starty_ = 0;
            
            for kk_ in 0..3{
                let kk = 2-kk_;
                for jj in 0..ylen{
                    if self.cells[xlen-1][jj].scores[kk] > maxscore{
                        maxscore = self.cells[xlen-1][jj].scores[kk];
                        startx_ = xlen as i64-1;
                        starty_ = jj as i64;
                        max_place = kk;
                    }
                }
            }
            for kk_ in 0..3{
                let kk = 2-kk_;
                for ii in 0..xlen{
                    if self.cells[ii][ylen-1].scores[kk] > maxscore{
                        maxscore = self.cells[ii][ylen-1].scores[kk];
                        startx_ = ii as i64;
                        starty_ = ylen as i64-1;
                        max_place = kk;
                    }
                }
            }
            */
        }
        if startx_ < 0{
            return (vec![],vec![],0.0);
        }
        let mut currentx:usize = startx_ as usize;
        let mut currenty:usize = starty_ as usize;
        let mut current_direc = max_place;
        loop{
            let prev_direc = self.cells[currentx][currenty].prev[current_direc];
            if self.alignment_type == ALIGN_LOCAL 
            && self.cells[currentx][currenty].scores[current_direc] == 0.0_f32{
                break;
            }
            //println!("direc:{} x:{} y:{}",current_direc,currentx,currenty);
            if current_direc == CELL_MATCH{
                ret1.push((currentx-1) as i64);
                ret2.push((currenty-1) as i64);
                currentx -= 1;
                currenty -= 1;
            }else if  current_direc == CELL_GAPINX{
                ret1.push(-1);
                ret2.push((currenty-1) as i64);
                currenty -= 1;
            }else if  current_direc == CELL_GAPINY{
                if currentx == 0{
                    panic!("{}",currenty);
                }
                ret1.push((currentx-1) as i64);
                ret2.push(-1);
                currentx -= 1;
            }else{
                panic!("???");
            }
            if currentx == 0 && currenty == 0{
                break;
            }
            current_direc = prev_direc;
            
        }
        ret1.reverse();
        ret2.reverse();
        return (ret1,ret2,maxscore);
    }
    
    pub fn align(&mut self,s1:&SeqData,s2:&SeqData,retain_all:bool)
    ->(Vec<String>,Vec<String>,f32){
        return self.align_partial(s1,s2,retain_all,None,false);
    }
    pub fn align_partial(&mut self,s1:&SeqData,s2:&SeqData,retain_all:bool,partial_region:Option<(usize,usize)>,score_only:bool)
    ->(Vec<String>,Vec<String>,f32){
        
        let (seq1,seq2):(Vec<usize>,Vec<usize>) = match partial_region{
            Some(x)=>{
                (self.scoring_matrix.seq_to_index(s1,Some(x.0)),self.scoring_matrix.seq_to_index(s2,Some(x.1)))
            },
            None => {
                (self.scoring_matrix.seq_to_index(s1,None),self.scoring_matrix.seq_to_index(s2,None))
            }
        };
        
        self.fill_matrix(&seq1,&seq2,partial_region);
        /*
        //スコア表を表示する部分 ===============
        let xlen = self.seqlen_a +1;
        let ylen = self.seqlen_b +1;
        
        for jj in 0..ylen{
            for ii in 0..xlen{
                print!("{} ",self.cells[ii][jj].scores[CELL_MATCH]);
            }
            println!("");
        }
        println!("");
        
        for jj in 0..ylen{
            for ii in 0..xlen{
                print!("{} ",self.cells[ii][jj].scores[CELL_GAPINX]);
            }
            println!("");
        }
        println!("");
        
        for jj in 0..ylen{
            for ii in 0..xlen{
                print!("{} ",self.cells[ii][jj].scores[CELL_GAPINY]);
            }
            println!("");
        }
        println!("");

        //スコア表を表示する部分終わり ===============
        */

        if score_only{
            let mut maxscore;
            if self.alignment_type == ALIGN_LOCAL{
                let xlen = self.seqlen_a +1;
                let ylen = self.seqlen_b +1;
                maxscore = 0.0;
                for ii in 0..xlen{
                    for jj in 0..ylen{
                        if self.cells[ii][jj].scores[CELL_MATCH] > maxscore{
                            maxscore = self.cells[ii][jj].scores[CELL_MATCH];
                        }
                    }
                }
            }else if self.alignment_type == ALIGN_GLOBAL || self.alignment_type == ALIGN_GLOCAL{
                let xlen = self.seqlen_a +1;
                let ylen = self.seqlen_b +1;
                let startx_ = xlen as i64-1;
                let starty_ = ylen as i64-1;
                let slen = self.cells[startx_ as usize][starty_ as usize].scores.len();
                maxscore = self.cells[startx_ as usize][starty_ as usize].scores[0];
                for ii in 1..slen{
                    if maxscore < self.cells[startx_ as usize][starty_ as usize].scores[ii]{
                        maxscore = self.cells[startx_ as usize][starty_ as usize].scores[ii];
                    }
                }
            }else{
                panic!();
            }
            return (vec![],vec![],maxscore);
        }

        let res = self.backtrack();
        let mut ret1:Vec<String> =vec![];
        let mut ret2:Vec<String> =vec![];
        
        let mut start1_:i64 = -1;
        let mut start2_:i64 = -1;
        let mut end1_:i64 = -1;
        let mut end2_:i64 = -1;
        
        //アラインした領域の要素を文字に直す
        for ii in res.0.iter(){
            if *ii > -1{
                if start1_ < 0{
                    start1_ = *ii;
                }
                ret1.push(s1.seq[*ii as usize].clone());
                end1_ = *ii;
            }else{
                ret1.push("-".to_string());
            }
        }
        
        for ii in res.1.iter(){
            if *ii > -1{
                if start2_ < 0{
                    start2_ = *ii;
                }
                ret2.push(s2.seq[*ii as usize].clone());
                end2_ = *ii;
            }else{
                ret2.push("-".to_string());
            }
        }
        if self.alignment_type == ALIGN_LOCAL && !retain_all{
            return (ret1,ret2,res.2);
        }
        
        if self.alignment_type != ALIGN_LOCAL && !retain_all{
            eprintln!("The glocal or global mode will retain all letters.");
        }
        
        let mut rret1:Vec<String> = vec![];
        let mut rret2:Vec<String> = vec![];
        if start1_ < 0 || start2_ < 0{
            //全く並ばなかった
            for ii in 0..s1.seq.len(){
                rret1.push(s1.seq[ii].clone());
                rret2.push("-".to_string());
            }    
            
            for ii in 0..s2.seq.len(){
                rret1.push("-".to_string());
                rret2.push(s2.seq[ii].clone());
            }  
            return (rret1,rret2,res.2);
        }
        //LOCAL でも全長を表示する
        let start1:usize = start1_ as usize;
        let start2:usize = start2_ as usize;
        let end1:usize = end1_ as usize;
        let end2:usize = end2_ as usize;
        for ii in 0..start1{
            rret1.push(s1.seq[ii].clone());
            rret2.push("-".to_string());
        }
        for ii in 0..start2{
            rret1.push("-".to_string());
            rret2.push(s2.seq[ii].clone());
        }
        
        rret1.append(&mut ret1);
        rret2.append(&mut ret2);
        
        for ii in (end1+1)..s1.seq.len(){
            rret1.push(s1.seq[ii].clone());
            rret2.push("-".to_string());
        }
        for ii in (end2+1)..s2.seq.len(){
            rret1.push("-".to_string());
            rret2.push(s2.seq[ii].clone());
        }
        return (rret1,rret2,res.2);
    }
}



#[derive(Clone)]
pub struct SWCell{
    pub scores:Vec<f32>,
    pub prev:Vec<usize>
}



impl SWCell{
    pub fn new()->SWCell{
        return SWCell{scores:vec![0.0;3],prev:vec![0_usize,0_usize,0_usize]};
    }
    pub fn set(&mut self,index:usize,prevpath:usize,score:f32){
        self.prev[index] = prevpath;
        self.scores[index] = score;
    }
}

pub trait ScoringMatrix{
    fn get_score(&self,a:usize,b:usize)->f32;
    fn get_score_str(&self,a:&str,b:&str)->f32;
    fn seq_to_index(&self,ss:&SeqData,partial_region:Option<usize>)->Vec<usize>;
    fn set_score(&mut self,a:usize,b:usize,s:f32);
    fn prepare(&mut self,a:&SeqData,b:&SeqData);
}

#[allow(dead_code)]
pub struct PositionSpecificMatrix{
    pub scores:Vec<f32>,
    pub a_length:usize,
    pub b_length:usize
}

impl ScoringMatrix for PositionSpecificMatrix{
    fn get_score_str(&self,_a:&str,_b:&str)->f32{
        panic!("not implemented");
    }
    fn get_score(&self,a:usize,b:usize)->f32{
        return self.scores[a+b*self.a_length];
    }
    fn seq_to_index(&self,ss:&SeqData,partial_region:Option<usize>)->Vec<usize>{
        if let Some(x) = partial_region{
        return (0..x).into_iter().collect();
        }else{
            return (0..ss.seq.len()).into_iter().collect();
        }
    }
    fn set_score(&mut self,a:usize,b:usize,s:f32){
        self.scores[a+b*self.a_length] = s;
    }
    fn prepare(&mut self,a:&SeqData,b:&SeqData){
        self.a_length = a.seq.len();
        self.b_length = b.seq.len();
        if self.a_length*self.b_length > self.scores.len(){
            self.scores = vec![0.0;self.a_length*self.b_length];
        }
    }
}

impl PositionSpecificMatrix{
    pub fn new()->PositionSpecificMatrix{
        return PositionSpecificMatrix{
            scores:vec![],
            a_length:0,
            b_length:0
        };
    }
}

#[allow(dead_code)]
pub struct SubstitutionMatrix{
    pub string_to_index:HashMap<String,usize> ,
    pub index_to_string:Vec<String>,
    pub scores:Vec<Vec<f32>>,
    
}
impl ScoringMatrix for SubstitutionMatrix{
    fn get_score_str(&self,a:&str,b:&str)->f32{
        
        if !self.string_to_index.contains_key(a){
            panic!("{} was not found in scoring matrix!",a);
        }
        if !self.string_to_index.contains_key(b){
            panic!("{} was not found in scoring matrix!",b);
        }

        return self.scores[*self.string_to_index.get(a).unwrap()
        ][*self.string_to_index.get(b).unwrap()
        ];
    }
    fn get_score(&self,a:usize,b:usize)->f32{
        return self.scores[a][b];
    }
    fn seq_to_index(&self,ss:&SeqData,partial_region:Option<usize>)->Vec<usize>{
        if let Some(x) = partial_region{
            let mut ret:Vec<usize> = vec![];
            for xx in 0..x{
                ret.push(self.get_string_index(&ss.seq[xx]));
            }
            return ret;
        }else{ 
            return ss.seq.iter().map(|m| self.get_string_index(m)).collect();
        }
    }
    fn set_score(&mut self,a:usize,b:usize,s:f32){
        self.scores[a][b] = s;
    }
    fn prepare(&mut self,_a:&SeqData,_b:&SeqData){
        //特に何もしない
    }
}

impl SubstitutionMatrix{
    pub fn get_string_index(&self,vv:&str)->usize{
        if self.string_to_index.contains_key(vv){
            return *(self.string_to_index.get(vv).unwrap());
        }else{
            if self.string_to_index.contains_key("X"){
                return *(self.string_to_index.get("X").unwrap());
            }else{
                panic!("unknown letter {}. please set X to allow scoring for undefined letter pair.");
            }
        }
    }
    //同アルファベットの場合 matchscore, 別アルファベットの場合 mismatchscore となるマトリクスを返す
    pub fn get_mat_matrix(matchscore:f32,mismatchscore:f32)->SubstitutionMatrix{
        let letters:Vec<String> = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".chars().map(|m|{m.to_string()}).collect();
        let mut string_to_index:HashMap<String,usize> = HashMap::new();
        let mut index_to_string:Vec<String> = vec![];
        let mut scores:Vec<Vec<f32>> = vec![vec![mismatchscore;letters.len()];letters.len()];
        for (ii,ll) in letters.into_iter().enumerate(){
            string_to_index.insert(ll.clone(),ii);
            index_to_string.push(ll);
            scores[ii][ii] = matchscore;
        }
        let ret:SubstitutionMatrix=SubstitutionMatrix{string_to_index:string_to_index,
        index_to_string:index_to_string,
        scores:scores};
        return ret;
    }

    pub fn get_blosum62_matrix()->SubstitutionMatrix{
        let mut lin:Vec<String> = Vec::new();
        //https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
        //#  Matrix made by matblas from blosum62.iij
        //#  * column uses minimum score
        //#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
        //#  Blocks Database = /data/blocks_5.0/blocks.dat
        //#  Cluster Percentage: >= 62
        //#  Entropy =   0.6979, Expected =  -0.5209
        lin.push("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *".to_string());
        lin.push("A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4".to_string());
        lin.push("R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4".to_string());
        lin.push("N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4".to_string());
        lin.push("D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4".to_string());
        lin.push("C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4".to_string());
        lin.push("Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4".to_string());
        lin.push("E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4".to_string());
        lin.push("G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4".to_string());
        lin.push("H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4".to_string());
        lin.push("I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4".to_string());
        lin.push("L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4".to_string());
        lin.push("K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4".to_string());
        lin.push("M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4".to_string());
        lin.push("F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4".to_string());
        lin.push("P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4".to_string());
        lin.push("S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4".to_string());
        lin.push("T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4".to_string());
        lin.push("W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4".to_string());
        lin.push("Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4".to_string());
        lin.push("V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4".to_string());
        lin.push("B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4".to_string());
        lin.push("Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4".to_string());
        lin.push("X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4".to_string());
        lin.push("* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1".to_string());

        return SubstitutionMatrix::load_from_lines(lin);
    }

    pub fn load_from_lines(lines:Vec<String>)-> SubstitutionMatrix{
        let mut scores_hm:HashMap<String,f32> = HashMap::new();
        let mut string_to_index:HashMap<String,usize>  = HashMap::new();
        let mut index_to_string:Vec<String>  = vec![];
        let mut column_to_c:HashMap<usize,String> = HashMap::new();//何番目のカラムが何という文字か
        let mut lincount:i64 = -1;
        for (_ii,line) in lines.into_iter().enumerate() {
            let bs = line.trim();
            let ptt:Vec<String> = bs.split_whitespace().map(|m| m.to_string()).collect();
            let lcc:char = ptt[0].chars().next().unwrap();
            if lcc == '#'{
                continue;
            }
            lincount+=1;
            if lincount == 0{
                for (jj,pp) in ptt.into_iter().enumerate(){
                    if string_to_index.contains_key(&pp){
                       panic!("{} was already found.",pp);
                    }
                    let csiz = string_to_index.keys().len();
                    string_to_index.insert(pp.clone(),csiz);
                    index_to_string.push(pp.clone());
                    column_to_c.insert(jj,pp.clone());
                }
            }else{
                if !string_to_index.contains_key(&ptt[0]){
                    panic!("{} was not found in the row name.",&ptt[0]);
                }

                for ll in 1..ptt.len(){
                    let ss = match ptt[ll].parse::<f32>(){
                        Ok(x)=>{
                            x
                        },
                        _=>{
                            eprintln!("{} can not be parsed! zero was assigned",ptt[ll]);
                            0.0
                        }
                    };
                    scores_hm.insert(format!("{} {}",lcc,column_to_c.get(&(ll-1)).unwrap()),ss);
                }
            }
        }
        let ilen = string_to_index.keys().len();
        let mut scores:Vec<Vec<f32>> =  vec![vec![0.0;ilen];ilen];
        for ii in 0..ilen{
            for jj in 0..ilen{
                let strcode = format!("{} {}",column_to_c.get(&ii).unwrap(),column_to_c.get(&jj).unwrap());
                if !scores_hm.contains_key(&strcode){
                    panic!("score about {} is not defined.",strcode);
                }
                scores[ii][jj] = scores_hm.get(&strcode).unwrap().clone();
            }
        }
        assert_eq!(index_to_string.len(),string_to_index.len());
        let ret:SubstitutionMatrix=SubstitutionMatrix{string_to_index:string_to_index,
        index_to_string:index_to_string,
        scores:scores};
        return ret;
    }
}

#[allow(dead_code)]
#[derive(Clone)]
pub struct SeqData{
    pub name:String,
    pub desc:String,
    pub seq:Vec<String>,
}
impl SeqData{
    pub fn new()->SeqData{
        return SeqData{name:"".to_string(),desc:"".to_string(),seq:Vec::new()};
    }
    pub fn create(n:String,d:String,s:String)->SeqData{
        return SeqData{name:n,desc:d,seq:SeqData::line_to_seq(s,true)};
    }
    pub fn line_to_seq(s:String,retainws:bool)->Vec<String>{
        //1 文字 1 アミノ酸を意味する文字列が与えられた場合
        if retainws {
            //pdb の ss.txt のパース目的
            let char_vec: Vec<char> = s.chars().collect();
            return char_vec.into_iter().filter(|m| m != &'\r' &&  m != &'\n' ).map(|m| m.to_string()).collect();
        }else{
            let char_vec: Vec<char> = s.chars().collect();
            return char_vec.into_iter().filter(|m| !m.is_whitespace() ).map(|m| m.to_string()).collect();
        }
    }

    /*
    pub fn filt_ws(&mut self){
        let mut vv:Vec<String> = Vec::new();
        vv.append(&mut self.seq);
        self.seq =  vv.into_iter().filter(|m| !m.is_whitespace() ).collect();
    }
    */
    pub fn load_fasta(filename:&str,retainws:bool)-> Vec<SeqData>{
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut ret:Vec<SeqData> = Vec::new();
        let mut seqbuff:Vec<String> = Vec::new();
        let mut currentname:String = "".to_string();
        let mut currentdesc:String = "".to_string();
        
        for (_lcount,line) in reader.lines().enumerate() {
            let line = line.unwrap();
            let rres = line.find(">");
            match rres{
                Some(x)=>{
                    //名前も配列もないものは無視される
                    if seqbuff.len() > 0 || currentname.len() > 0{
                        let mut vv:Vec<String> = Vec::new();
                        vv.append(&mut seqbuff);//中身だけ移動
                        ret.push(SeqData{name:currentname.clone(),desc:currentdesc.clone(),seq:vv});
                    }
                    if x > 0{
                        eprintln!("> was found at {}. This line was used as header anyway.",x);
                    }
                    let line = line.trim();
                    let mut nameflag = true;
                    let mut namebuff:Vec<char> = Vec::new();
                    let mut descbuff:Vec<char> = Vec::new();
                    for (sii,sss) in line.chars().enumerate(){
                        if nameflag{
                            if sii == 0 && sss == '>'{
                                continue;
                            }
                            if sss.is_whitespace(){
                                if namebuff.len()>0{
                                    nameflag = false;
                                }
                                continue;
                            }
                            namebuff.push(sss);
                        }else{
                            descbuff.push(sss);
                        }
                    }
                    currentname = namebuff.iter().collect();
                    currentdesc = descbuff.iter().collect();
                },
                _=>{
                    seqbuff.append(&mut SeqData::line_to_seq(line,retainws));
                }
            }

        }
        
        if currentname.len() > 0 ||  seqbuff.len() > 0{
            ret.push(SeqData{name:currentname,desc:currentdesc,seq:seqbuff});
        }

        return ret;
    }
    
}






#[test]
fn sw_scoringmatrixtest(){
    let mut vlin:Vec<String> = Vec::new();
    //ncbi によれば default は Match +1 Mismatch -3 らしい
    //https://www.ebi.ac.uk/Tools/sss/ncbiblast/help/index-vectors.html
    /*
    Match/mismatch_scores
    (Nucleotide searches) The match score is the bonus to the alignment score when matching the same base. The mismatch is the penalty when failing to match.
    */
    vlin.push("  A T G C \n".to_string());
    vlin.push("A 1 -3 -3 -3 ".to_string());
    vlin.push("T -3  1 -3 -3 ".to_string());
    vlin.push("G -3  -3  1 -3 ".to_string());
    vlin.push("C -3  -3  -3  1".to_string());
    let sm = SubstitutionMatrix::load_from_lines(vlin);
    assert_eq!(sm.get_score_str("A","A"),1.0);
    assert_eq!(sm.get_score_str("T","T"),1.0);
    assert_eq!(sm.get_score_str("G","G"),1.0);
    assert_eq!(sm.get_score_str("C","C"),1.0);
    assert_eq!(sm.get_score_str("A","C"),-3.0);
    assert_eq!(sm.get_score_str("T","G"),-3.0);
    assert_eq!(sm.get_score_str("G","A"),-3.0);
    assert_eq!(sm.get_score_str("C","T"),-3.0);
}
#[test]
fn sw_aligntest(){
    let sm = SubstitutionMatrix::get_mat_matrix(5.0_f32,-4.0_f32);
    //">NC_024511.2:3768-3838 Drosophila melanogaster mitochondrion, complete genome"
    let seq1_ = "CATTAGATGACTGAAAGCAAGTACTGGTCTCTTAAACCATTTAATAGTAAATTAGCACTTACTTCTAATGA".to_string();
    //">NC_024462.2:196885816-196885888 Zea mays cultivar B73 chromosome 4, B73 RefGen_v4, whole genome shotgun sequence"
    let seq2_ = "ACTTCTCTAGCTCAGTTGGTAGAGCGCAAGGCTTTTAACCTTGTGGTCGTGGGTTCAAACCCCATGATGGGCA".to_string();

    let mut sw = SequenceAlignment::new(Box::new(sm),10.0,0.5,ALIGN_LOCAL);

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    //ebi の EMBOSS で得た結果
    assert_eq!(res.2,54.5);
    assert_eq!(r1,"CAT-------TAGATGACT-----GAAAGCAAG----------TACTGGTC------TCTTAAACCATTTAATAGTAAATTAGCACTTACTTCTAATGA");
    assert_eq!(r2,"---ACTTCTCTAGCTCAGTTGGTAGAGCGCAAGGCTTTTAACCTTGTGGTCGTGGGTTC--AAACCCCATGATGG-------GCA--------------");

    let sm = SubstitutionMatrix::get_mat_matrix(5.0_f32,-4.0_f32);
    let mut sw = SequenceAlignment::new(Box::new(sm),10.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    
    assert_eq!(res.2,51.5);
    assert_eq!(r1,"----CATTAGATGACT-----GAAAGCAAG----------TACTGGTC------TCTTAAACCATTTAATAGTAAATTAGCACTTACTTCTAATGA");
    assert_eq!(r2,"ACTTCTCTAGCTCAGTTGGTAGAGCGCAAGGCTTTTAACCTTGTGGTCGTGGGTTC--AAACCCCATGATGG-------GCA--------------");

    let sm = SubstitutionMatrix::get_mat_matrix(5.0_f32,-4.0_f32);
    let mut sw = SequenceAlignment::new(Box::new(sm),10.0,0.5,ALIGN_GLOBAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,24.0);
    assert_eq!(r1,"CATT---------AGATGACTGAAAGCAAGTACTGGTCTCTTAAACCATTTAATAGTAAATTAGCACTTACTTC-TAATG---A");
    assert_eq!(r2,"ACTTCTCTAGCTCAGTTGGTAGAGCGCAAG-GCT------TTTAACCTTGTGGTCGTGGGTT----CAAACCCCATGATGGGCA");
}


#[test]
fn sw_aligntest2(){
    let seq1_ = "MASSWKLMLFLSVTMCLSEYSKSLPGLSTSYAALLRIKKSSSSSLFGSKTRPRYSSPSLGTLSASSPSWLGAAQNYYSPINLYHSSDAFKQDESVDYGPVFVQEPDDIIFPTDSDEKKVALNCEVRGNPVPSYRWLRNGTEIDLESDYRYSLIDGTFIISNPSEAKDSGHYQCLATNTVGSILSREATLQFAYLGNFSGRTRSAVSVREGQGVVLMCSPPPHSPEIIYSWVFNEFPSFVAEDSRRFISQETGNLYISKVQTSDVGSYICLVKNTVTNARVLSPPTPLTLRNDGVMGEYEPKIEVHFPFTVTAAKGTTVKMECFALGNPVPTITWMKVNGYIPSKARLRKSQAVLEIPNVQLDDAGIYECRAENSRGKNSFRGQLQVYTYPHWVEKLNDTQLDSGSPLRWECKATGKPRPTYRWLKNGVPLSPQSRVEMVNGVLMIHNVNQSDAGMYQCLAENKYGAIYASAELKILASAPTFALNQLKKTIIVTKDQEVVIECKPQGSPKPTISWKKGDRAVRENKRIAILPDGSLRILNASKSDEGKYVCRGENVFGSAEIIASLSVKEPTRIELTPKRTELTVGESIVLNCKAIHDASLDVTFYWTLKGQPIDFEEEGGHFESIRAQASSADLMIRNILLMHAGRYGCRVQTTADSVSDEAELLVRGPPGPPGIVIVEEITESTATLSWSPAADNHSPISSYNLQARSPFSLGWQTVKTVPEIITGDMESAMAVDLNPWVEYEFRVVATNPIGTGDPSTPSRMIRTNEAVPKTAPTNVSGRSGRRHELVIAWEPVSEEFQNGEGFGYIVAFRPNGTRGWKEKMVTSSEASKFIYRDESVPPLTPFEVKVGVYNNKGDGPFSQIVVICSAEGEPSAAPTDVKATSVSVSEILVAWKHIKESLGRPQGFEVGYWKDMEQEDTAETVKTRGNESFVILTGLEGNTLYHFTVRAYNGAGYGPPSSEVSATTKKSPPSQAPSNLRWEQQGSQVSLGWEPVIPLANESEVVGYKVFYRQEGHSNSQVIETQKLQAVVPLPDAGVYIIEVRAYSEGGDGTASSQIRVPSYSGGKITSAQSTLHSLSTSSSSVTLLLALMIPSTSW".to_string();
    let seq2_ = "MTIRLLCYVGFYFLGAGLMEADIYQTPRYLVIGTGKKITLECSQTMGHDKMYWYQQDPGMELHLIHYSYGVNSTEKGDLSSESTVSRIRTEHFPLTLESARPSHTSQYLCASSE".to_string();

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let sm = SubstitutionMatrix::get_blosum62_matrix();
    let mut sw = SequenceAlignment::new(Box::new(sm),10.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,23.0);
    assert_eq!(r1,"MASSWKLMLFLSVTMCLSEYSKSLPGLSTSYAALLRIKKSSSSSLFGSKTRPRYSSPSLGTLSASSPSWLGAAQNYYSPINLYHSSDAFKQDESVDYGPVFVQEPDDIIFPTDSDEKKVALNCEVRGNPVPSYRWLRNGTEIDLESDYRYSLIDGTFIISNPSEAKDSGHYQCLATNTVGSILSREATLQ-FAYLG-NFSG---------RTRSAVSVREGQGVVLMCS--------------PPPHSPEIIYSWVFN--EFPSFVAEDSRRFISQETGNLYISKVQTSDVGSYICLVKNTVTNARVLSPPTPLTLRNDGVMGEYEPKIEVHFPFTVTAAKGTTVKMECFALGNPVPTITWMKVNGYIPSKARLRKSQAVLEIPNVQLDDAGIYECRAENSRGKNSFRGQLQVYTYPHWVEKLNDTQLDSGSPLRWECKATGKPRPTYRWLKNGVPLSPQSRVEMVNGVLMIHNVNQSDAGMYQCLAENKYGAIYASAELKILASAPTFALNQLKKTIIVTKDQEVVIECKPQGSPKPTISWKKGDRAVRENKRIAILPDGSLRILNASKSDEGKYVCRGENVFGSAEIIASLSVKEPTRIELTPKRTELTVGESIVLNCKAIHDASLDVTFYWTLKGQPIDFEEEGGHFESIRAQASSADLMIRNILLMHAGRYGCRVQTTADSVSDEAELLVRGPPGPPGIVIVEEITESTATLSWSPAADNHSPISSYNLQARSPFSLGWQTVKTVPEIITGDMESAMAVDLNPWVEYEFRVVATNPIGTGDPSTPSRMIRTNEAVPKTAPTNVSGRSGRRHELVIAWEPVSEEFQNGEGFGYIVAFRPNGTRGWKEKMVTSSEASKFIYRDESVPPLTPFEVKVGVYNNKGDGPFSQIVVICSAEGEPSAAPTDVKATSVSVSEILVAWKHIKESLGRPQGFEVGYWKDMEQEDTAETVKTRGNESFVILTGLEGNTLYHFTVRAYNGAGYGPPSSEVSATTKKSPPSQAPSNLRWEQQGSQVSLGWEPVIPLANESEVVGYKVFYRQEGHSNSQVIETQKLQAVVPLPDAGVYIIEVRAYSEGGDGTASSQIRVPSYSGGKITSAQSTLHSLSTSSSSVTLLLALMIPSTSW");
    assert_eq!(r2,"------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MTIRLLCYVGFYFLGAGLMEADIYQTPRYLVIGTGKKITLECSQTMGHDKMYWYQQDPGMELHLIHYSYGVNSTEKGDLSSESTVSRIRTEHFPLTLESARPSHTSQYLCASSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
    

    
    
    let seq1_ = "RRHRSEDCGGGPRSLSRGLPCKKAATEGSSEKTVLDSKPSVPTTSEGGPELELQIPELPLDSNEFWVHEGCILWANGIYLVCGRLYGLQEALEIAREMKCSHCQEAGATLGCYNKGCSFRYHYPCAIDADCLLHEENFSVRCPKHKPPLPCPLPPLQNKTAKGSLSTEQSERG".to_string();
    let seq2_ = "MKLAFLFLGPMALLLLAGYGCVLGASSGNLRTFVGCAVREFTFLAKKPGCRGLRITTDACWGRCETWEKPILEPPYIEAHHRVCTYNETKQVTVKLPNCAPGVDPFYTYPVAIRCDCGACSTATTECETI".to_string();
    
    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    //let sm = ScoringMatrix::get_blosum62_matrix();
    //let mut sw = SequenceAlignment::new(sm,10.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,22.0);
    assert_eq!(r1,"RRHRSEDCGG--------GPRS--LSRGLPCKKAATEGS---------SEKTVLDSKPSVPTTSEGGPELELQIPELPLDSNEFWVHEGCILWANGIYLVCGRLYGLQEALEIAREMKCSHCQEAGATLGCYN--KGCSFRYHYPCAIDADCLLHEENFSVRCPKHKPPLPCPLPPLQNKTAKGSLSTEQSERG--");
    assert_eq!(r2,"----------MKLAFLFLGPMALLLLAGYGCVLGASSGNLRTFVGCAVREFTFLAKKPG----CRG----------LRITTDACWGR--CETWEKPI---------LEPPYIEAHHRVCTYNETKQVTVKLPNCAPGVDPFYTYPVAIRCDC-------------------------------GACSTATTECETI");
    


    let seq1_ = "MQKIMHISVLLSPVLWGLIFGVSSNSIQIGGLFPRGADQEYSAFRVGMVQFSTSEFRLTPHIDNLEVANSFAVTNAFCSQFSRGVYAIFGFYDKKSVNTITSFCGTLHVSFITPSFPTDGTHPFVIQMRPDLKGALLSLIEYYQWDKFAYLYDSDRGLSTLQAVLDSAAEKKWQVTAINVGNINNDKKDEMYRSLFQDLELKKERRVILDCERDKVNDIVDQVITIGKHVKGYHYIIANLGFTDGDLLKIQFGGANVSGFQIVDYDDSLVSKFIERWSTLEEKEYPGAHTTTIKYTSALTYDAVQVMTEAFRNLRKQRIEISRRGNAGDCLANPAVPWGQGVEIERALKQVQVEGLSGNIKFDQNGKRINYTINIMELKTNGPRKIGYWSEVDKMVVTLTELPSGNDTSGLENKTVVVTTILESPYVMMKKNHEMLEGNERYEGYCVDLAAEIAKHCGFKYKLTIVGDGKYGARDADTKIWNGMVGELVYGKADIAIAPLTITLVREEVIDFSKPFMSLGISIMIKKPQKSKPGVFSFLDPLAYEIWMCIVFAYIGVSVVLFLVSRFSPYEWHTEEFEDGRETQSSESTNEFGIFNSLWFSLGAFMQQGCDISPRSLSGRIVGGVWWFFTLIIISSYTANLAAFLTVERMVSPIESAEDLSKQTEIAYGTLDSGSTKEFFRRSKIAVFDKMWTYMRSAEPSVFVRTTAEGVARVRKSKGKYAYLLESTMNEYIEQRKPCDTMKVGGNLDSKGYGIATPKGSSLRNAVNLAVLKLNEQGLLDKLKNKWWYDKGECGSGGGDSKEKTSALSLSNVAGVFYILVGGLGLAMLVALIEFCYKSRAEAKRMKVAKNAQNINPSSSQNSQNFATYKEGYNVYGIESVKI".to_string();
    let seq2_ = "MAPGSRTSLLLAFALLCLPWLQEAGAVQTVPLSRLFDHAMLQAHRAHQLAIDTYQEFEETYIPKDQKYSFLHDSQTSFCFSDSIPTPSNMEETQQKSNLELLRISLLLIESWLEPVRFLRSMFANNLVYDTSDSDDYHLLKDLEEGIQTLMGRLEDGSRRTGQILKQTYSKFDTNSHNHDALLKNYGLLYCFRKDMDKVETFLRMVQCRSVEGSCGF".to_string();
    
    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    //let sm = ScoringMatrix::get_blosum62_matrix();
    //let mut sw = SequenceAlignment::new(sm,10.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(r1,"MQKIMHISVLLSPVLWGLIFGVSSNSIQ---IGGLFPRGADQEYSAFRVGMVQFSTSEFRLTPHIDNLEVANSFAVTNAFCSQFSRGVYAIFGFYDKKSVNTITSFCGTLHVSFITPSFPTDGTHPFVIQMRPDLKGALLSLIEYYQW----------------------DKFAYLYDSDRGLSTLQAVLDSAAEKKWQV-----TAINVGNINNDKKDEMYRSLFQDLELKKERRVILDCER---DKVNDIVDQVITIGKHVKGYHYIIANLGFTDGDLLKIQFGGANVSGFQIVDYDDSLVSKFIERWSTLEEKEYPGAHTTTIKYTSALTYDAVQVMTEAFRNLRKQRIEISRRGNAGDCLANPAVPWGQGVEIERALKQVQVEGLSGNIKFDQNGKRINYTINIMELKTNGPRKIGYWSEVDKMVVTLTELPSGNDTSGLENKTVVVTTILESPYVMMKKNHEMLEGNERYEGYCVDLAAEIAKHCGFKYKLTIVGDGKYGARDADTKIWNGMVGELVYGKADIAIAPLTITLVREEVIDFSKPFMSLGISIMIKKPQKSKPGVFSFLDPLAYEIWMCIVFAYIGVSVVLFLVSRFSPYEWHTEEFEDGRETQSSESTNEFGIFNSLWFSLGAFMQQGCDISPRSLSGRIVGGVWWFFTLIIISSYTANLAAFLTVERMVSPIESAEDLSKQTEIAYGTLDSGSTKEFFRRSKIAVFDKMWTYMRSAEPSVFVRTTAEGVARVRKSKGKYAYLLESTMNEYIEQRKPCDTMKVGGNLDSKGYGIATPKGSSLRNAVNLAVLKLNEQGLLDKLKNKWWYDKGECGSGGGDSKEKTSALSLSNVAGVFYILVGGLGLAMLVALIEFCYKSRAEAKRMKVAKNAQNINPSSSQNSQNFATYKEGYNVYGIESVKI");
    assert_eq!(r2,"MAPGSRTSLLLAFALLCLPWLQEAGAVQTVPLSRLFDHAMLQAHRAHQLAIDTYQEFEETYIPK----DQKYSF-------------------LHDSQ-----TSFC-------FSDSIPTP-SNMEETQQKSNLELLRISLLLIESWLEPVRFLRSMFANNLVYDTSDSDDYHLLKDLEEGIQTLMGRLEDGSRRTGQILKQTYSKFDTNSHNHDALLKNYGLLY--------------CFRKDMDKVETFLRMVQC--RSVEG------SCGF-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
    
    
    
    let seq1_ = "MALQDVCKWQSPDTQGPSPHLPRAGGWAVPRGCDPQTFLQIHGPRLAHGTTTLAFRFRHGVIAAADTRSSCGSYVACPASCKVIPVHQHLLGTTSGTSADCATWYRVLQRELRLRELREGQLPSVASAAKLLSAMMSQYRGLDLCVATALCGWDRSGPELFYVYSDGTRLQGDIFSVGSGSPYAYGVLDRGYRYDMSTQEAYALARCAVAHATHRDAYSGGSVDLFHVRESGWEHVSRSDACVLYVELQKLLEPEPEEDASHAHPEPATAHRAAEDRELSVGPGEVTPGDSRMPAGTETV".to_string();
    let seq2_ = "MAVAPSFNMTNPQPAIEGGISEVEIISQQVDEETKSIAPVQLVNFAYRDLPLAAVDLSTAGSQLLSNLDEDYQREGSNWLKPCCGKRAAVWQVFLLSASLNSFLVACVILVVILLTLELLIDIKLLQFSSAFQFAGVIHWISLVILSVFFSETVLRIVVLGIWDYIENKIEVFDGAVIILSLAPMVASTVANGPRSPWDAISLIIMLRIWRVKRVIDAYVLPVKLEMEMVIQQYEKAKVIQDEQLERLTQICQEQGFEIRQLRAHLAQQDLDLAAEREAALQAPHVLSQPRSRFKVLEAGTWDEETAAESVVEELQPSQEATMKDDMNSYISQYYNGPSSDSGVPEPAVCMVTTAAIDIHQPNISSDLFSLDMPLKLGGNGTSATSESASRSSVTRAQSDSSQTLGSSMDCSTAREEPSSEPGPSPPPLPSQQQVEEATVQDLLSSLSEDPCPSQKALDPAPLARPSPAGSAQTSPELEHRVSLFNQKNQEGFTVFQIRPVIHFQPTVPMLEDKFRSLESKEQKLHRVPEA".to_string();

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let sm = SubstitutionMatrix::get_blosum62_matrix();
    let mut sw = SequenceAlignment::new(Box::new(sm),10.0,0.5,ALIGN_GLOBAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(r1,"MALQDVCKWQSPDTQGPSPHLPRAGGWAVPRGCDPQTFLQIHGPRLAHGTTTLA------FRFRHGVIAAADTRSSCGSYVAC---------------PASCKVIPVHQHLLGTTSGTS--ADCATWYRVLQRELRLRELREGQLPSVASAAK-----------LLSAMMSQYRGLDLCVATALCG-WD--RSGPELFYVYSDGTRLQGDIFSV----------GSGSPY-------------------AYGVLDRGYRYDMSTQEAYALAR-------------CAV---------AHATHRDAYSGGSVDL--------------------FHVRESG-WEHVSRSDACVLYVELQKLLEPEPEED-------------ASHAHPEPA--------------------------------------------TAHRAAEDRELSVG---------------PG-------------EVT-------------PGDSRM---------PAGTE---------------------TV------------------------------------");
    assert_eq!(r2,"MAV-------APSFNMTNPQPAIEGGISE---------VEIISQQVDEETKSIAPVQLVNFAYRDLPLAAVDL-STAGSQLLSNLDEDYQREGSNWLKPCCGKRAAVWQVFLLSASLNSFLVACVILVVIL---LTLELLIDIKLLQFSSAFQFAGVIHWISLVILSVFFS-----ETVLRIVVLGIWDYIENKIEVF----DGAVI---ILSLAPMVASTVANGPRSPWDAISLIIMLRIWRVKRVIDAY-VLPVKLEMEMVIQQ-YEKAKVIQDEQLERLTQICQEQGFEIRQLRAHLAQQD------LDLAAEREAALQAPHVLSQPRSRFKVLEAGTWDEETAAESVV--EELQPSQEATMKDDMNSYISQYYNGPSSDSGVPEPAVCMVTTAAIDIHQPNISSDLFSLDMPLKLGGNGTSATSESASRSSVTRAQSDSSQTLGSSMDCSTAREEPSSEPGPSPPPLPSQQQVEEATVQDLLSSLSEDPCPSQKALDPAPLARPSPAGSAQTSPELEHRVSLFNQKNQEGFTVFQIRPVIHFQPTVPMLEDKFRSLESKEQKLHRVPEA");
    
    
    let seq1_ = "MQKIMHISVLLSPVLWGLIFGVSSNSIQIGGLFPRGADQEYSAFRVGMVQFSTSEFRLTPHIDNLEVANSFAVTNAFCSQFSRGVYAIFGFYDKKSVNTITSFCGTLHVSFITPSFPTDGTHPFVIQMRPDLKGALLSLIEYYQWDKFAYLYDSDRGLSTLQAVLDSAAEKKWQVTAINVGNINNDKKDEMYRSLFQDLELKKERRVILDCERDKVNDIVDQVITIGKHVKGYHYIIANLGFTDGDLLKIQFGGANVSGFQIVDYDDSLVSKFIERWSTLEEKEYPGAHTTTIKYTSALTYDAVQVMTEAFRNLRKQRIEISRRGNAGDCLANPAVPWGQGVEIERALKQVQVEGLSGNIKFDQNGKRINYTINIMELKTNGPRKIGYWSEVDKMVVTLTELPSGNDTSGLENKTVVVTTILESPYVMMKKNHEMLEGNERYEGYCVDLAAEIAKHCGFKYKLTIVGDGKYGARDADTKIWNGMVGELVYGKADIAIAPLTITLVREEVIDFSKPFMSLGISIMIKKPQKSKPGVFSFLDPLAYEIWMCIVFAYIGVSVVLFLVSRFSPYEWHTEEFEDGRETQSSESTNEFGIFNSLWFSLGAFMQQGCDISPRSLSGRIVGGVWWFFTLIIISSYTANLAAFLTVERMVSPIESAEDLSKQTEIAYGTLDSGSTKEFFRRSKIAVFDKMWTYMRSAEPSVFVRTTAEGVARVRKSKGKYAYLLESTMNEYIEQRKPCDTMKVGGNLDSKGYGIATPKGSSLRNAVNLAVLKLNEQGLLDKLKNKWWYDKGECGSGGGDSKEKTSALSLSNVAGVFYILVGGLGLAMLVALIEFCYKSRAEAKRMKVAKNAQNINPSSSQNSQNFATYKEGYNVYGIESVKI".to_string();
    let seq2_ = "MAPGSRTSLLLAFALLCLPWLQEAGAVQTVPLSRLFDHAMLQAHRAHQLAIDTYQEFEETYIPKDQKYSFLHDSQTSFCFSDSIPTPSNMEETQQKSNLELLRISLLLIESWLEPVRFLRSMFANNLVYDTSDSDDYHLLKDLEEGIQTLMGRLEDGSRRTGQILKQTYSKFDTNSHNHDALLKNYGLLYCFRKDMDKVETFLRMVQCRSVEGSCGF".to_string();
    
    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    //let sm = ScoringMatrix::get_blosum62_matrix();
    //let mut sw = SequenceAlignment::new(sm,10.0,0.5,ALIGN_GLOBAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(r1,"MQKIMHISVLLSPVLWGLIFGVSSNSIQ---IGGLFPRGADQEYSAFRVGMVQFSTSEFRLTPHIDNLEVANSFAVTNAFCSQFSRGVYAIFGFYDKKSVNTITSFCGTLHVSFITPSFPTDGTHPFVIQMRPDLKGALLSLIEYYQWDKFAYLYDSDRGLSTLQAVLDSAAEKKWQVTAINVGNINNDKKDEMYRSLFQDLELKKERRVILDCERDKVNDIVDQVITIGKHVKGYHYIIANLGFTDGDLLKIQFGGANVSGFQIVDYDDSLVSKFIERWSTLEEKEYPGAHTTTIKYTSALTYDAVQVMTEAFRNLRKQRIEISRRGNAGDCLANPAVPWGQGVEIERALKQVQVEGLSGNIKFDQNGKRINYTINIMELKTNGPRKIGYWSEVDKMVVTLTELPSGNDTSGLENKTVVVTTILESPYVMMKKNHEMLEGNERYEGYCVDLAAEIAKHCGFKYKLTIVGDGKYGARDADTKIWNGMVGELVYGKADIAIAPLTITLVREEVIDFSKPFMSLGISIMIKKPQKSKPGVFSFLDPLAYEIWMCIVFAYIGVSVVLFLVSRFSPYEWHTEEFEDGRETQSSESTNEFGIFNSLWFSLGAFMQQGCDISPRSLSGRIVGGVWWFFTLIIISSYTANLAAFLTVERMVSPIESAEDLSKQTEIAYGTLDSGSTKEFFRRSKIAVFDKMWTYMRSAEPSVFVRTTAEGVARVRKSKGKYAYLLESTMNEYIEQRKPCDTMKVGGNLDSKGYGIATPKGSSLRNAVNLAVLKLNEQGLLDKLKNKWWYDKGECGSGGGDSKEKTSALSLSNVAGVFYILVGGLGLAMLVALIEFCYKSRAEAKRMKVAKNAQNINPSSSQNSQNFATYKEGYNVYGIESVKI");
    assert_eq!(r2,"MAPGSRTSLLLAFALLCLPWLQEAGAVQTVPLSRLFDHAMLQAHRAHQLAIDTYQEFEETYIPK----DQKYSF-------------------LHDSQ-----TSFC-------FSDSIPTP-SNMEETQQKSNLELLRISLLLIESW---------------LEPV-------------------------RFLRSMF-------------------ANNLV---------------------------------------YDTSDSDD---------------------------------YHLLKDLEEGIQTLMGRLEDGSRR-------------TGQ------ILKQTY-------SKFDTN----------------------------------------------------------------SHNHDALLKN-------------------------------YG---------------LLY----------------------------------------------------------C----------------------------------------------------------------------------------------------------------------------------------FRKD----MDKVETFLRMVQ---------------------------------------CRSV------------------------------------------------EGSCG---------------------------------------------------------------------------------------F");
    

    
    let seq1_ = "MISPDPRPSPGLARWAESYEAKCERRQEIRESRRCRPNVTTCRQVGKTLRIQQREQLQRARLQQFFRRRNLELEEKGKAQHPQAREQGPSRRPGQVTVLKEPLSCARRISSPREQVTGTSSEVFPAQHPPPSGICRDLSDHLSSQAGGLPPQDTPIKKPPKHHRGTQTKAEGPTIKNDASQQTNYGVAVLDKEIIQLSDYLKEALQRELVLKQKMVILQDLLSTLIQASDSSWKGQLNEDKLKGKLRSLENQLYTCTQKYSPWGMKKVLLEMEDQKNSYEQKAKESLQKVLEEKMNAEQQLQSTQRSLALAEQKCEEWRSQYEALKEDWRTLGTQHRELESQLHVLQSKLQGADSRDLQMNQALRFLENEHQQLQAKIECLQGDRDLCSLDTQDLQDQLKRSEAEKLTLVTRVQQLQGLLQNQSLQLQEQEKLLTKKDQALPVWSPKSFPNEVEPEGTGKEKDWDLRDQLQKKTLQLQAKEKECRELHSELDNLSDEYLSCLRKLQHCREELNQSQQLPPRRQCGRWLPVLMVVIAAALAVFLANKDNLMI".to_string();
    let seq2_ = "MNPTETKAIPVSQQMEGPHLPNKKKHKKQAVKTEPEKKSQSTKLSVVHEKKSQEGKPKEHTEPKSLPKQASDTGSNDAHNKKAVSRSAEQQPSEKSTEPKTKPQDMISAGGESVAGITAISGKPGDKKKEKKSLTPAVPVESKPDKPSGKSGMDAALDDLIDTLGGPEETEEENTTYTGPEVSDPMSSTYIEELGKREVTIPPKYRELLAKKEGITGPPADSSKPIGPDDAIDALSSDFTCGSPTAAGKKTEKEESTEVLKAQSAGTVRSAAPPQEKKRKVEKDTMSDQALEALSASLGTRQAEPELDLRSIKEVDEAKAKEEKLEKCGEDDETIPSEYRLKPATDKDGKPLLPEPEEKPKPRSESELIDELSEDFDRSECKEKPSKPTEKTEESKAAAPAPVSEAVCRTSMCSIQSAPPEPATLKGTVPDDAVEALADSLGKKEADPEDGKPVMDKVKEKAKEEDREKLGEKEETIPPDYRLEEVKDKDGKPLLPKESKEQLPPMSEDFLLDALSEDFSGPQNASSLKFEDAKLAAAISEVVSQTPASTTQAGAPPRDTSQSDKDLDDALDKLSDSLGQRQPDPDENKPMEDKVKEKAKAEHRDKLGERDDTIPPEYRHLLDDNGQDKPVKPPTKKSEDSKKPADDQDPIDALSGDLDSCPSTTETSQNTAKDKCKKAASSSKAPKNGGKAKDSAKTTEETSKPKDD".to_string();
    
        
    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    //let sm = ScoringMatrix::get_blosum62_matrix();
    //let mut sw = SequenceAlignment::new(sm,10.0,0.5,ALIGN_GLOBAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(r1,"MISPDPRPSPGLARWAESY--EAKCERRQEIRESRRCRPNVTTCRQVGKTLRIQQREQLQRARLQQFFRRRNL----------------------ELEEKGKAQHPQAREQ--------------GPSRRPGQVTVLKEPLSCARRISSPREQVTGTSS---------------EVFPAQHPPPSGICRDLSDHLSS----------------------QAGGL--PPQDT--PI----------------------KKPPKHHRGTQTKAEGP-TIKNDASQQTNYGVAVLDKEIIQLSDYLKEALQRELVLKQKMVILQDLLSTLIQASDSSWKGQLNEDKLK-----GKLRSLENQLYTCTQKYSPWGMKKVLLEMEDQKNSYEQKAKESLQKVLEE------KMNAEQQLQSTQRSLALA-----EQKCEEWRSQYEALKEDWRTL-GTQHRELESQLHVLQSKLQG--ADSRDLQMNQALRFLENEHQQLQAKIECLQGDRDLCSLDTQDLQDQLKRSEAEKLTLVTRVQQLQGLL------QNQSLQLQEQEKLLTKKDQALPVWSPKSFPNEVEPEGTGKEKDWDLRDQLQK---------------KTLQLQAKEKECRELHSEL----DNLSDEYLSCLRKLQHCREELNQSQQLPP---------------------------------------RRQCGRWLPVLMVVIAAALAVFLAN----KDNLMI----------");
    assert_eq!(r2,"MNPTETKAIPVSQQMEGPHLPNKKKHKKQAVKTEPEKKSQST-------KLSVVHEKKSQEGKPKEHTEPKSLPKQASDTGSNDAHNKKAVSRSAEQQPSEKSTEPKTKPQDMISAGGESVAGITAISGKPGDKKKEKKSLTPAVPVESKPDKPSGKSGMDAALDDLIDTLGGPEETEEENTTYTG--PEVSDPMSSTYIEELGKREVTIPPKYRELLAKKEGITGPPADSSKPIGPDDAIDALSSDFTCGSPTAAGKKTEKEESTEVLKAQSAGTVRSAAPPQEKKRKVEKD----TMSDQALEALSASLGTRQAEPEL-DLRS--IKEVD---EAKAKEEKLEKCGEDDETIPSEYRLKPATDKDG----KPLLPEPEEKPK--PRSESELIDELSEDFDRSECKEKPSKPTEKTEESKAAAPAPVSEAVCRTSMCSIQSAPPEPATLKGTVPDDAVEALADSLGKKEADPEDGKPVMDKVKEKAKEEDREKLGEKEETIPPDYRLEEVKDKDGKPLLPKESKEQLPPMSEDFLLDALSEDFSGPQNASSLKFEDAKLAAAISEVVS-QTPASTTQAGAPPRDTSQSDKDLDDALDKLSDSLGQRQPDPDENKPMEDKVKEKAKAEHRDKLGERDDTIPPEY-------RHLLDDNGQDKPVKPPTKKSEDSKKPADDQDPIDALSGDLDSCPSTTETSQNTAKDKCKK---------AASSSKAPKNGGKAKDSAKTTEETSKPKDD");

}

#[test]
fn psm_test(){
    


    let seq1_ = "MASSWKLMLFLSVTMCLSEYSKSLPGLSTSYAALLRIKKSSSSSLFGSKTRPRYSSPSLGTLSASSPSWLGAAQNYYSPINLYHSSDAFKQDESVDYGPVFVQEPDDIIFPTDSDEKKVALNCEVRGNPVPSYRWLRNGTEIDLESDYRYSLIDGTFIISNPSEAKDSGHYQCLATNTVGSILSREATLQFAYLGNFSGRTRSAVSVREGQGVVLMCSPPPHSPEIIYSWVFNEFPSFVAEDSRRFISQETGNLYISKVQTSDVGSYICLVKNTVTNARVLSPPTPLTLRNDGVMGEYEPKIEVHFPFTVTAAKGTTVKMECFALGNPVPTITWMKVNGYIPSKARLRKSQAVLEIPNVQLDDAGIYECRAENSRGKNSFRGQLQVYTYPHWVEKLNDTQLDSGSPLRWECKATGKPRPTYRWLKNGVPLSPQSRVEMVNGVLMIHNVNQSDAGMYQCLAENKYGAIYASAELKILASAPTFALNQLKKTIIVTKDQEVVIECKPQGSPKPTISWKKGDRAVRENKRIAILPDGSLRILNASKSDEGKYVCRGENVFGSAEIIASLSVKEPTRIELTPKRTELTVGESIVLNCKAIHDASLDVTFYWTLKGQPIDFEEEGGHFESIRAQASSADLMIRNILLMHAGRYGCRVQTTADSVSDEAELLVRGPPGPPGIVIVEEITESTATLSWSPAADNHSPISSYNLQARSPFSLGWQTVKTVPEIITGDMESAMAVDLNPWVEYEFRVVATNPIGTGDPSTPSRMIRTNEAVPKTAPTNVSGRSGRRHELVIAWEPVSEEFQNGEGFGYIVAFRPNGTRGWKEKMVTSSEASKFIYRDESVPPLTPFEVKVGVYNNKGDGPFSQIVVICSAEGEPSAAPTDVKATSVSVSEILVAWKHIKESLGRPQGFEVGYWKDMEQEDTAETVKTRGNESFVILTGLEGNTLYHFTVRAYNGAGYGPPSSEVSATTKKSPPSQAPSNLRWEQQGSQVSLGWEPVIPLANESEVVGYKVFYRQEGHSNSQVIETQKLQAVVPLPDAGVYIIEVRAYSEGGDGTASSQIRVPSYSGGKITSAQSTLHSLSTSSSSVTLLLALMIPSTSW".to_string();
    let seq2_ = "MTIRLLCYVGFYFLGAGLMEADIYQTPRYLVIGTGKKITLECSQTMGHDKMYWYQQDPGMELHLIHYSYGVNSTEKGDLSSESTVSRIRTEHFPLTLESARPSHTSQYLCASSE".to_string();


    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let pm = PositionSpecificMatrix::new();
    let sm = SubstitutionMatrix::get_blosum62_matrix();
    let svec1:Vec<usize> = sm.seq_to_index(&seq1,None) ;
    let svec2:Vec<usize> = sm.seq_to_index(&seq2,None) ;
    let mut sw = SequenceAlignment::new(Box::new(pm),10.0,0.5,ALIGN_GLOCAL);
    sw.prepare(&seq1,&seq2);
    for ii in 0..svec1.len(){
        for jj in 0..svec2.len(){
            sw.scoring_matrix.set_score(ii,jj, sm.get_score(svec1[ii],svec2[jj]));
        }
    }
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,23.0);
    assert_eq!(r1,"MASSWKLMLFLSVTMCLSEYSKSLPGLSTSYAALLRIKKSSSSSLFGSKTRPRYSSPSLGTLSASSPSWLGAAQNYYSPINLYHSSDAFKQDESVDYGPVFVQEPDDIIFPTDSDEKKVALNCEVRGNPVPSYRWLRNGTEIDLESDYRYSLIDGTFIISNPSEAKDSGHYQCLATNTVGSILSREATLQ-FAYLG-NFSG---------RTRSAVSVREGQGVVLMCS--------------PPPHSPEIIYSWVFN--EFPSFVAEDSRRFISQETGNLYISKVQTSDVGSYICLVKNTVTNARVLSPPTPLTLRNDGVMGEYEPKIEVHFPFTVTAAKGTTVKMECFALGNPVPTITWMKVNGYIPSKARLRKSQAVLEIPNVQLDDAGIYECRAENSRGKNSFRGQLQVYTYPHWVEKLNDTQLDSGSPLRWECKATGKPRPTYRWLKNGVPLSPQSRVEMVNGVLMIHNVNQSDAGMYQCLAENKYGAIYASAELKILASAPTFALNQLKKTIIVTKDQEVVIECKPQGSPKPTISWKKGDRAVRENKRIAILPDGSLRILNASKSDEGKYVCRGENVFGSAEIIASLSVKEPTRIELTPKRTELTVGESIVLNCKAIHDASLDVTFYWTLKGQPIDFEEEGGHFESIRAQASSADLMIRNILLMHAGRYGCRVQTTADSVSDEAELLVRGPPGPPGIVIVEEITESTATLSWSPAADNHSPISSYNLQARSPFSLGWQTVKTVPEIITGDMESAMAVDLNPWVEYEFRVVATNPIGTGDPSTPSRMIRTNEAVPKTAPTNVSGRSGRRHELVIAWEPVSEEFQNGEGFGYIVAFRPNGTRGWKEKMVTSSEASKFIYRDESVPPLTPFEVKVGVYNNKGDGPFSQIVVICSAEGEPSAAPTDVKATSVSVSEILVAWKHIKESLGRPQGFEVGYWKDMEQEDTAETVKTRGNESFVILTGLEGNTLYHFTVRAYNGAGYGPPSSEVSATTKKSPPSQAPSNLRWEQQGSQVSLGWEPVIPLANESEVVGYKVFYRQEGHSNSQVIETQKLQAVVPLPDAGVYIIEVRAYSEGGDGTASSQIRVPSYSGGKITSAQSTLHSLSTSSSSVTLLLALMIPSTSW");
    assert_eq!(r2,"------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MTIRLLCYVGFYFLGAGLMEADIYQTPRYLVIGTGKKITLECSQTMGHDKMYWYQQDPGMELHLIHYSYGVNSTEKGDLSSESTVSRIRTEHFPLTLESARPSHTSQYLCASSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
    println!("{}",r2);

}


#[test]
fn sw_fastaloadtest(){
    let fas = SeqData::load_fasta("test/test1.fas",false);
    assert_eq!(fas[0].name,"seqA".to_string());
    assert_eq!(fas[0].desc,"".to_string());
    assert_eq!(fas[0].seq,vec!["A";12]);
    
    assert_eq!(fas[1].name,"seqB".to_string());
    assert_eq!(fas[1].desc,"b desu".to_string());
    assert_eq!(fas[1].seq,vec!["B";12]);
    
    assert_eq!(fas[2].name,"seqC".to_string());
    assert_eq!(fas[2].desc,"c desu".to_string());
    assert_eq!(fas[2].seq,vec!["C";12]);

    assert_eq!(fas[3].name,"seqD".to_string());
    assert_eq!(fas[3].desc,"d desu".to_string());
    let emptyv:Vec<String> = vec![];
    assert_eq!(fas[3].seq,emptyv);
    

    assert_eq!(fas[4].name,"seqE".to_string());
    assert_eq!(fas[4].desc,"e desu".to_string());
    assert_eq!(fas[4].seq,vec!["E";12]);
    
    assert_eq!(fas[5].name,"F".to_string());
    assert_eq!(fas[5].desc,"".to_string());
    assert_eq!(fas[5].seq,vec!["F";12]);
    
    assert_eq!(fas[6].name,"G".to_string());
    assert_eq!(fas[6].desc,"".to_string());
    assert_eq!(fas[6].seq,vec!["G";12]);
    
    assert_eq!(fas[7].name,"H".to_string());
    assert_eq!(fas[7].desc,"h".to_string());
    assert_eq!(fas[7].seq,vec!["H";12]);
    

    assert_eq!(fas[8].name,"I".to_string());
    assert_eq!(fas[8].desc,"".to_string());
    assert_eq!(fas[8].seq,emptyv);
}


#[test]
fn sw_aligntest3(){
    //Glocal のスコアは違うかもしれない
    let seq1_ = "AAAAASSSSSS".to_string();
    let seq2_ = "NNNNNSSSSSS".to_string();

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let sm = SubstitutionMatrix::get_blosum62_matrix();
    let mut sw = SequenceAlignment::new(Box::new(sm),8.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,14.0);
    //assert_eq!(r1,"AAAAASSSSSS");
    //assert_eq!(r2,"NNNNNSSSSSS");

    let seq1_ = "AAAAASSSSSS".to_string();
    let seq2_ = "NNNNNSSSSSS".to_string();

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let sm = SubstitutionMatrix::get_blosum62_matrix();
    let mut sw = SequenceAlignment::new(Box::new(sm),7.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,15.0);
    assert_eq!(r1,"AAAAA-----SSSSSS");
    assert_eq!(r2,"-----NNNNNSSSSSS");

    let seq1_ = "SSSSSSAAAAA".to_string();
    let seq2_ = "SSSSSSNNNNN".to_string();

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let sm = SubstitutionMatrix::get_blosum62_matrix();
    let mut sw = SequenceAlignment::new(Box::new(sm),7.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,15.0);
    assert_eq!(r1,"SSSSSSAAAAA-----");
    assert_eq!(r2,"SSSSSS-----NNNNN");
    
    
    let seq1_ = "MKVSEAALSLLVLILIITSASRSQPKVPEWVNTPSTCCLKYYEKVLPRRLVVGYRKALNCHLPAIIFVTKRNREVCTNPNDDWVQEYIKDPNLPLLPTRNLSTVKIITAKNGQPQLLNSQ".to_string();
    let seq2_ = "MGNITADNSSMSCTIDHTIHQTLAPVVYVTVLVVGFPANCLSLYFGYLQIKARNELGVYLCNLTVADLFYICSLPFWLQYVLQHDNWSHGDLSCQVCGILLYENIYISVGFLCCISVDRYLAVAHPFRFHQFRTLKAAVGVSVVIWAKELLTSIYFLMHEEVIEDENQHRVCFEHYPIQAWQRAINYYRFLVGFLFPICLLLASYQGILRAVRRSHGTQKSRKDQIQRLVLSTVVIFLACFLPYHVLLLVRSVWEASCDFAKGVFNAYHFSLLLTSFNCVADPVLYCFVSETTHRDLARLRGACLAFLTCSRTGRAREAYPLGAPEASGKSGAQGEEPELLTKLHPAFQTPNSPGSGGFPTGRLA".to_string();

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let sm = SubstitutionMatrix::get_blosum62_matrix();
    let mut sw = SequenceAlignment::new(Box::new(sm),10.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,10.5);
    assert_eq!(r1,"-----------------------------------------------------------------------------------------------------------------------------------MKVSEAALSLLVLIL---IITS----------ASRSQPKVPEWVNTPSTCCLKYYEKVLPRRLVVGYRKALNCHLPAIIFV---------------TKRNREVCTNPNDDWVQEYIKDPNLPLLPTRNLSTVKIITAKNGQPQLLNSQ------------------------------------------------------------------------------------------------------------------------------");
    assert_eq!(r2,"MGNITADNSSMSCTIDHTIHQTLAPVVYVTVLVVGFPANCLSLYFGYLQIKARNELGVYLCNLTVADLFYICSLPFWLQYVLQHDNWSHGDLSCQVCGILLYENIYISVGFLCCISVDRYLAVAHPFRFHQFRTLKAAVGVSVVIWAKELLTSIYFLMHEEVIEDENQHRV----------CFEHYPIQAWQRAINYYRFLVGFLFPICLLLASYQGILRAVRRSHGTQKSRK-------DQIQRLV------------LSTVVIFLA-----------CFLPYHVLLLVRSVWEASCDFAKGVFNAYHFSLLLTSFNCVADPVLYCFVSETTHRDLARLRGACLAFLTCSRTGRAREAYPLGAPEASGKSGAQGEEPELLTKLHPAFQTPNSPGSGGFPTGRLA");

    let seq1_ = "MESLRGYTHSDIGYRSLAVGEDIEEVNDEKLTVTSLMARGGEDEENTRSKPEYGTEAENNVGTEGSVPSDDQDREGGGGHEPEQQQEEPPLTKPEQQQEEPPLLELKQEQEEPPQTTVEGPQPAEGPQTAEGPQPPERKRRRRTAFTQFQLQELENFFDESQYPDVVARERLAARLNLTEDRVQVWFQNRRAKWKRNQRVLMLRNTATADLAHPLDMFLGGAYYAAPALDPALCVHLVPQLPRPPVLPVPPMPPRPPMVPMPPRPPIAPMPPMAPVPPGSRMAPVPPGPRMAPVPPWPPMAPVPPWPPMAPVPTGPPMAPVPPGPPMARVPPGPPMARVPPGPPMAPLPPGPPMAPLPPGPPMAPLPPGPPMAPLPPRSHVPHTGLAPVHITWAPVINSYYACPFF".to_string();
    let seq2_ = "MPNVLLPPKESNLFKRILKCYEQKQYKNGLKFCKMILSNPKFAEHGETLAMKGLTLNCLGKKEEAYEFVRKGLRNDVKSHVCWHVYGLLQRSDKKYDEAIKCYRNALKLDKDNLQILRDLSLLQIQMRDLEGYRETRYQLLQLRPTQRASWIGYAIAYHLLKDYDMALKLLEEFRQTQQVPPNKIDYEYSELILYQNQVMREADLLQESLEHIEMYEKQICDKLLVEEIKGEILLKLGRLKEASEVFKNLIDRNAENWCYYEGLEKALQISTLEERLQIYEEISKQHPKAITPRRLPLTLVPGERFRELMDKFLRVNFSKGCPPLFTTLKSLYYNTEKVSIIQELVTNYEASLKTCDFFSPYENGEKEPPTTLLWVQYFLAQHFDKLGQYSLALDYINAAIASTPTLIELFYMKAKIYKHIGNLKEAAKWMDEAQSLDTADRFINSKCAKYMLRANMIKEAEEMCSKFTREGTSAMENLNEMQCMWFQTECISAYQRLGRYGDALKKCHEVERHFFEITDDQFDFHTYCMRKMTLRAYVDLLRLEDILRRHAFYFKAARSAIEIYLKLYDNPLTNESKQQEINSENLSAKELKKMLSKQRRAQKKAKLEEERKHAERERQQKNQKKKRDEEEEEASGLKEELIPEKLERVENPLEEAVKFLIPLKNLVADNIDTHLLAFEIYFRKGKFLLMLQSVKRAFAINSNNPWLHECLIRFSKSVSNHSNLPDIVSKVLSQEMQKIFVKKDLESFNEDFLKRNATSLQHLLSGAKMMYFLDKSRQEKAIAIATRLDETIKDKDVKTLIKVSEALLDGSFGNCSSQYEEYRMACHNLLPFTSAFLPAVNEVDNPNVALNHTANYDVLANEI".to_string();

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let sm = SubstitutionMatrix::get_blosum62_matrix();
    let mut sw = SequenceAlignment::new(Box::new(sm),10.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,4.0);
    assert_eq!(r1,"MESLRGYTHSDIGYRSLAVGEDIEEVNDEKLTVTSLMARGGEDEENTRSKPEYGTEAENNVGTEGSVPSDDQDREGGGGHEPEQQQEEPPLTKPEQQQEEPPLLELKQEQEEPPQTTVEGPQPAEGPQTAEGPQPPERKRRRRTAFTQFQLQELENFFDESQYPDVVARERLAARLNLTEDRVQVWFQNRRAKWKRNQRVLMLRNTATADLAHPLDMFLGGAYYAAPALDPALCVHLVPQLPRPPVLPVPPMPPRPPMVPMPPRPPIAPMPPMAPVPPGSRMAPVPPGPRMAPVPPWPPMAPVPPWPPMAPVPTGPPMAPVPPGPPMARVPPGPPMARVPPGPPMAPLPPGPPMAPLPPGPPMAPLPPGPPMAPLPPRSHVPHTGLAPVHITWAPVINSYYACPFF----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
    assert_eq!(r2,"--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MPNVLLPPKESNLFKRILKCYEQKQYKNGLKFCKMILSNPKFAEHGETLAMKGLTLNCLGKKEEAYEFVRKGLRNDVKSHVCWHVYGLLQRSDKKYDEAIKCYRNALKLDKDNLQILRDLSLLQIQMRDLEGYRETRYQLLQLRPTQRASWIGYAIAYHLLKDYDMALKLLEEFRQTQQVPPNKIDYEYSELILYQNQVMREADLLQESLEHIEMYEKQICDKLLVEEIKGEILLKLGRLKEASEVFKNLIDRNAENWCYYEGLEKALQISTLEERLQIYEEISKQHPKAITPRRLPLTLVPGERFRELMDKFLRVNFSKGCPPLFTTLKSLYYNTEKVSIIQELVTNYEASLKTCDFFSPYENGEKEPPTTLLWVQYFLAQHFDKLGQYSLALDYINAAIASTPTLIELFYMKAKIYKHIGNLKEAAKWMDEAQSLDTADRFINSKCAKYMLRANMIKEAEEMCSKFTREGTSAMENLNEMQCMWFQTECISAYQRLGRYGDALKKCHEVERHFFEITDDQFDFHTYCMRKMTLRAYVDLLRLEDILRRHAFYFKAARSAIEIYLKLYDNPLTNESKQQEINSENLSAKELKKMLSKQRRAQKKAKLEEERKHAERERQQKNQKKKRDEEEEEASGLKEELIPEKLERVENPLEEAVKFLIPLKNLVADNIDTHLLAFEIYFRKGKFLLMLQSVKRAFAINSNNPWLHECLIRFSKSVSNHSNLPDIVSKVLSQEMQKIFVKKDLESFNEDFLKRNATSLQHLLSGAKMMYFLDKSRQEKAIAIATRLDETIKDKDVKTLIKVSEALLDGSFGNCSSQYEEYRMACHNLLPFTSAFLPAVNEVDNPNVALNHTANYDVLANEI");

    let seq1_ = "MSNATLLTAFILTGLPHAPGLDAPLFGIFLVVYVLTVLGNLLILLVIRVDSHLHTPMYYFLTNLSFIDMWFSTVTVPKMLMTLVSPSGRTISFHSCVAQLYFFHFLGSTECFLYTVMSYDRYLAISYPLRYTNMMTGRSCALLATGTWLSGSLHSAVQTILTFHLPYCGPNQIQHYFCDAPPILKLACADTSANEMVIFVNIGLVASGCFVLIVLSYVSIVCSILRIRTSEGRHRAFQTCASHCIVVLCFFGPGLFIYLRPGSRDALHGVVAVFYTTLTPLFNPVVYTLRNKEVKKALLKLKNGSVFAQGE".to_string();
    let seq2_ = "MMNNTDFLMLNNPWNKLCLVSMDFCFPLDFVSNLFWIFASKFIIVTGQIKADFKRTSWEAKAEGSLEPGRLKLQLASIVPLYSSLVTAGPASKIIILKRTSLPTVSPSNERAYLLPVSFTDLAHVFYLSYFSINAKSNSFSLDIIIALGIPHNTQAHFNH".to_string();

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    let sm = SubstitutionMatrix::get_blosum62_matrix();
    let mut sw = SequenceAlignment::new(Box::new(sm),10.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(res.2,19.5);
    assert_eq!(r1,"---------------------------------------------------------------------------------------------------------------------------------------MSNATLLTAFILTGLP--------HAPGLDAPLFGIFLVVYVLTVLGNLLILLVIRVDSHLHTPMYYFLTNLSFIDMWFSTVTVPKMLMTLVSPSGRTISFHSCVAQLYFFHFLGSTECFLYTVMSYDRYLAISYPLRYTNMMTGRSCALLATGTWLSGSLHSAVQTILTFHLPYCGPNQIQHYFCDAPPILKLACADTSANEMVIFVNIGLVASGCFVLIVLSYVSIVCSILRIRTSEGRHRAFQTCASHCIVVLCFFGPGLFIYLRPGSRDALHGVVAVFYTTLTPLFNPVVYTLRNKEVKKALLKLKNGSVFAQGE");
    assert_eq!(r2,"MMNNTDFLMLNNPWNKLCLVSMDFCFPLDFVSNLFWIFASKFIIVTGQIKADFKRTSWEAKAEGSLEPGRLKLQLASIVPLYSSLVTAGPASKIIILKRTSLPTVSPSNERAYLLPVSFTDLAHVFYLSYFSINAKSNSFSLDIIIALGIPHNTQAHFNH------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");

}

