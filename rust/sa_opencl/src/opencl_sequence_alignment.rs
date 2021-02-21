
use std::fs::File;
use std::collections::HashMap;
use std::io::{BufReader,BufRead};
use ocl::{Platform, Device, Context, Queue, Program,Buffer, Kernel};


// https://github.com/yamule/smithwaterman
// OpenCL Smith-Waterman & Needleman-Wunsch sequence alignment 
// using anti-diagonal parallelization, which is similar with the approach
// Wozniak, Andrzej. "Using video-oriented instructions to speed up sequence comparison." Bioinformatics 13.2 (1997): 145-150.
// Awan, Muaaz G., et al. "ADEPT: a domain independent sequence alignment strategy for gpu architectures." BMC bioinformatics 21.1 (2020): 1-29.
//
// According to 
// Awan, Muaaz G., et al. "ADEPT: a domain independent sequence alignment strategy for gpu architectures." BMC bioinformatics 21.1 (2020): 1-29.
// ,
// Edans, F. de O., et al. "CUDAlign 3.0: Parallel biological sequence comparison in large GPU clusters." 2014 14th IEEE/ACM International Symposium on Cluster, Cloud and Grid Computing. IEEE, 2014.
// also uses the approach but I was blocked by the paywall.
//
// License: Public Domain

const POS_TERMINAL:i32 = -999;
pub struct VecAndBuffer<T:ocl_core::OclPrm>{
    pub vector:Vec<T>,
    pub buffer:ocl::Buffer<T>,
}
impl<T:ocl_core::OclPrm> VecAndBuffer<T>{
    pub fn new(vvec:Vec<T>,qque:&ocl::Queue)->VecAndBuffer<T>{
        let buffer =  match 
        ocl::Buffer::<T>::builder()
        .queue(qque.clone())
        .flags(ocl::MemFlags::new().read_write().copy_host_ptr())
        .len(vvec.len())
        .copy_host_slice(&vvec.as_slice())
        .build(){
            Ok(x) => {x},
            Err(e) =>{
                panic!("{:?}",e);
            }
        };
        
        if let Err(x) = buffer.write(&vvec).enq(){
            panic!("{:?}",x);
        }
        return VecAndBuffer{vector:vvec,buffer:buffer};
    }

    pub fn read(&mut self){
        if let Err(x) = self.buffer.read(&mut self.vector).enq(){
            panic!("error when read from buff {:?}",x);
        }
    }
    pub fn write(&mut self){
        if let Err(x) = self.buffer.write(&self.vector).enq(){
            panic!("error when read from buff {:?}",x);
        }
    }

}



pub const CELL_MATCH:usize = 0;
pub const CELL_GAPINX:usize = 1;
pub const CELL_GAPINY:usize = 2;

pub const ALIGN_GLOBAL:usize = 0;
pub const ALIGN_GLOCAL:usize = 1;
pub const ALIGN_LOCAL:usize = 2;

pub struct OpenCLSequenceAlignment{
    pub scoring_matrix:Box<dyn ScoringMatrix>,
    pub scoring_matrix_vec:Buffer<f32>,
    pub e_penalty:f32,
    pub o_penalty:f32,
    pub start_epenal:f32,
    pub start_openal:f32,
    pub seq_length:VecAndBuffer<i32>,
    pub num_cols:usize,//== seqlen_a+1
    pub num_rows:usize,//== seqlen_b+1
    pub alignment_type:usize,
    pub kernel_fill:Kernel,
    //pub kernel_backtrack:Kernel,

    pub ocl_platform:Platform,
    pub ocl_device:Device,
    pub ocl_context:Context,
    pub ocl_program:Program,
    pub ocl_queue:Queue,
    pub seq_a:VecAndBuffer<i32>,
    pub seq_b:VecAndBuffer<i32>,
    pub aligned_seq_a:VecAndBuffer<i32>,
    pub aligned_seq_b:VecAndBuffer<i32>,
    pub max_position:VecAndBuffer<i32>,
    pub score_buff:VecAndBuffer<f32>,
    pub dp_matrix:VecAndBuffer<f32>,
    pub flag_matrix:VecAndBuffer<u8>,
    pub kernel_control:VecAndBuffer<u8>,

}

impl OpenCLSequenceAlignment{
    pub fn new(max_length:usize,ss:Box<dyn ScoringMatrix>,go:f32,ge:f32,alignment_type:usize)->OpenCLSequenceAlignment{
        let src_test:String = format!(r#"
        int pos_2d_to_1d_rc(int r,int c,int num_rows,int num_cols){{
            
            int k = max((c+r-num_rows+1),0);
            int l = max((c+r-num_cols+1),0);
            int ks = k*num_rows;
            int ls = l*num_cols;
            int pt = l*k;
            int ps = ks+ls-pt;
            int cc = c-k;
            int rr = r-l;
            int zt = (cc+rr)*(cc+rr+1)/2+rr;

            return ps+zt;
        }}
        int pos_2d_to_1d_rc_(int r,int c,int num_rows,int num_cols){{
            return num_cols*r+c;
        }}
        
        __kernel void prepare_matrix(
            __global float* dp_matrix
            ,__global char* flag_matrix // 00 match 01 gapincol 10 gapincol x3 0 uncalculated 1 calculated
            ,__global int* seqsize
            ,float start_openal
            ,float start_epenal) {{
            
            int cc = get_global_id(0);
            int num_rows = seqsize[1]+1;
            int num_cols = seqsize[0]+1;
            if(get_global_id(0) >= num_cols){{return;}}
            for(int rr = 0;rr < num_rows;rr++){{
                int ppos = pos_2d_to_1d_rc(rr,cc,num_rows,num_cols);
                int ppos_d = ppos*3;
                flag_matrix[ppos] = 0b0000000;
                dp_matrix[ppos_d] = -10000.0;
                dp_matrix[ppos_d+1] = -10000.0;
                dp_matrix[ppos_d+2] = -10000.0;
            }}
    
            for(int pp_ = 0;pp_ < 2;pp_++){{
                int tc = (pp_ == 0)?(num_cols):(1);
                int tr = (pp_ == 1)?(num_rows):(1);
                for(int rr = 0;rr < tr;rr++){{
                    float openal = start_openal;
                    float epenal = start_epenal;
                    int ppos = pos_2d_to_1d_rc(rr,cc,num_rows,num_cols);
                    int ppos_d = ppos*3;

                    if(cc == 0 && rr == 0){{
                        dp_matrix[ppos_d] = 0.0;
                        dp_matrix[ppos_d+1] = 0.0;
                        dp_matrix[ppos_d+2] = 0.0;
                        flag_matrix[ppos] = 0b0000001;
                        continue;
                    }}
                    if(cc == 0){{
                        float lscore = rr*epenal+(openal-epenal);
                        dp_matrix[ppos_d] = lscore+10.0*openal +10.0*openal -10000.0;
                        dp_matrix[ppos_d+1] = lscore;
                        dp_matrix[ppos_d+2] = lscore+10.0*openal +10.0*openal -10000.0;
                        flag_matrix[ppos] = 0b0101011;
                        continue;
                    }}else if(rr == 0){{
                        float lscore = cc*epenal+(openal-epenal);
                        dp_matrix[ppos_d] = lscore+10.0*openal +10.0*openal -10000.0 ;
                        dp_matrix[ppos_d+1] = lscore+10.0*openal +10.0*openal -10000.0;
                        dp_matrix[ppos_d+2] = lscore;
                        flag_matrix[ppos] = 0b1010101;
                        continue;
                    }}
                }}
            }}
        }}
        __global void backtrack(
            __global int* aligned_seq_a 
            ,__global int* aligned_seq_b
            ,__global int* max_position
            ,__global float* score_buff
            ,__global int* seqsize
            ,__global float* dp_matrix
            ,__global char* flag_matrix // 00 match 01 gapincol 10 gapincol x3 0 uncalculated 1 calculated
            ,int alignment_type){{
            int startx_ = -1;
            int starty_ = -1;
            int num_rows = seqsize[1]+1;
            int num_cols = seqsize[0]+1;
            float maxscore = 0.0;
            int max_place = 0;
            score_buff[0] = 0.0;
            if(alignment_type == {ALIGN_LOCAL}){{
                maxscore = 0.0;
                for(int cc = 1;cc < num_cols;cc++){{
                    int rr = max_position[cc];
                    int ppos = pos_2d_to_1d_rc(rr,cc,num_rows,num_cols);
                    if(dp_matrix[ppos*3] > maxscore){{
                        maxscore = dp_matrix[ppos*3];
                        startx_ = cc;
                        starty_ = rr;
                    }}
                
                }}
                max_place = {CELL_MATCH};
            }}else if(alignment_type == {ALIGN_GLOBAL} || alignment_type == {ALIGN_GLOCAL}){{
                
                startx_ = num_cols -1;
                starty_ = num_rows -1;
                
                int ppos = pos_2d_to_1d_rc(starty_,startx_,num_rows,num_cols);

                maxscore = dp_matrix[ppos*3];
                max_place = 0;
                for(int ii = 0;ii <3;ii++){{
                    if(maxscore < dp_matrix[ppos*3+ii]){{
                        maxscore = dp_matrix[ppos*3+ii];
                        max_place = ii;
                    }}
                }}
            }}else{{
                printf("!!!!!!!!!!!!!!!!!!!!!!!!!!! this alignment type is not expected!!!!!!!");
                return;
            }}
            if(startx_ < 0){{
                aligned_seq_a[0] = {POS_TERMINAL};
                aligned_seq_b[0] = {POS_TERMINAL};
                return;
            }}
            int currentx = startx_;
            int currenty = starty_;
            int current_direc = max_place;
            
            int buffposx = 0;
            int buffposy = 0;
            while(1){{
                int currentpos = pos_2d_to_1d_rc(currenty,currentx,num_rows,num_cols);
                //printf("%d,",currentpos);
                char prev_direc_ = flag_matrix[currentpos];
                char prev_direc = 0;
                /*デバッグ中
                if((prev_direc_ & 0b1) == 0 ){{
                    printf("????pos :%d direc: %d????? %d %d %d %d \n",currentpos,current_direc,currentx,currenty,buffposx,buffposy);
                    //return;
                }}else{{
                    printf("!!!pos :%d direc: %d????? %d %d %d %d \n",currentpos,current_direc,currentx,currenty,buffposx,buffposy);
                }}
                */
                if(current_direc == {CELL_MATCH}){{
                    prev_direc = (prev_direc_ >> 1) & 0b11;
                }}else if(current_direc == {CELL_GAPINX}){{
                    prev_direc = (prev_direc_ >> 3) & 0b11;
                }}else if(current_direc == {CELL_GAPINY}){{
                    prev_direc = (prev_direc_ >> 5) & 0b11;
                }}else{{
                    printf("????This direction is not expected pos :%d direc: %d????? %d %d",currentpos,current_direc,currentx,currenty);
                    return;
                }};
                
                if (alignment_type == {ALIGN_LOCAL} 
                && dp_matrix[currentpos*3+current_direc] <= 0.0){{
                    break;
                }}
                
                if(current_direc == {CELL_MATCH}){{
                    aligned_seq_a[buffposx] = currentx-1;
                    aligned_seq_b[buffposy] = currenty-1;
                    currentx -= 1;
                    currenty -= 1;
                }}else if(current_direc == {CELL_GAPINX}){{
                    aligned_seq_a[buffposx] = -1;
                    aligned_seq_b[buffposy] = currenty-1;
                    currenty -= 1;
                }}else if(current_direc == {CELL_GAPINY}){{
                    aligned_seq_a[buffposx] = currentx-1;
                    aligned_seq_b[buffposy] = -1;
                    currentx -= 1;
                }}else{{
                    printf("???");
                    return;
                }}
                buffposx++;
                buffposy++;
                
                if(currentx == 0 && currenty == 0){{
                    break;
                }}
                if(currentx < 0 || currenty < 0){{
                    printf("!!! pos :%d direc: %d????? %d %d",currentpos,current_direc,currentx,currenty);
                    break;
                }}
                current_direc = prev_direc;
            }}
            
            aligned_seq_a[buffposx] = {POS_TERMINAL};
            aligned_seq_b[buffposy] = {POS_TERMINAL};
            score_buff[0] = maxscore;
        }}
        __kernel void fill_matrix(
            __global int* seq_a 
            ,__global int* seq_b
            ,__global int* seqsize
            ,__global float* scoring_matrix
            ,int num_letter_types //文字が何種類あるか
            ,__global float* dp_matrix
            ,__global char* flag_matrix // 00 match 01 gapincol 10 gapincol x3 0 uncalculated 1 calculated
            ,float start_openal
            ,float start_epenal
            ,float openal
            ,float epenal
            ,__global int* aligned_seq_a 
            ,__global int* aligned_seq_b
            ,__global int* max_position
            ,__global float* score_buff
            ,int alignment_type
            ,__global char* kernel_control//0 prepare, 1 fill, 2 backtrack
            
        ) {{
            if(kernel_control[0] == 0){{
                //printf("%d started %d\n",get_global_id(0),get_local_size(1));
                prepare_matrix(
                dp_matrix
                ,flag_matrix
                ,seqsize
                ,start_openal
                ,start_epenal);
                //printf("%d blocked\n",get_global_id(0));
                //barrier(CLK_GLOBAL_MEM_FENCE);
                //printf("%d running\n",get_global_id(0));
            }}else if(kernel_control[0] == 1){{
                int num_rows = seqsize[1]+1;
                int num_cols = seqsize[0]+1;

                int col_id = get_global_id(0)+1;
                int maxpos = -1;
                float maxscore = -10000;
                if(col_id < num_cols){{
                    int prevcol = col_id-1;
                    int prevrow = 0;
                    int currentrow = 1;
                    
                    int prevpos_t = pos_2d_to_1d_rc(currentrow-1,col_id,num_rows,num_cols);
                    int prevpos_l = pos_2d_to_1d_rc(currentrow,col_id-1,num_rows,num_cols);
                    int prevpos_lt = pos_2d_to_1d_rc(currentrow-1,col_id-1,num_rows,num_cols);
                    int currentpos = pos_2d_to_1d_rc(currentrow,col_id,num_rows,num_cols);
                    maxscore = dp_matrix[currentpos];
                    maxpos = currentrow;

                    while(currentrow < num_rows){{
                        if(
                            (flag_matrix[prevpos_t] & 1) == 1
                            &&
                            (flag_matrix[prevpos_l] & 1) == 1
                            &&
                            (flag_matrix[prevpos_lt] & 1) == 1
                        ){{

                            int p3t = prevpos_t*3;
                            int p3l = prevpos_l*3;
                            int p3lt = prevpos_lt*3;
                            int cp3 = currentpos*3;

                            float mmscore = scoring_matrix[seq_a[col_id-1]*num_letter_types+seq_b[currentrow-1]];
                            
                            float matchscore = 0.0;
                            int matchindex = 1;
                            
                            float gapxscore = 0.0;
                            int gapxindex = 0;
                            
                            float gapyscore = 0.0;
                            int gapyindex = 0;
                            
                            //Cell にマッチした場合
                            if(dp_matrix[p3lt+{CELL_MATCH}]+mmscore
                            >= dp_matrix[p3lt+{CELL_GAPINX}]+mmscore){{
                                if(dp_matrix[p3lt+{CELL_MATCH}]+mmscore
                                >= dp_matrix[p3lt+{CELL_GAPINY}]+mmscore){{    
                                    matchscore = dp_matrix[p3lt+{CELL_MATCH}]+mmscore;
                                    matchindex = {CELL_MATCH};
                                }}else{{
                                    matchscore = dp_matrix[p3lt+{CELL_GAPINY}]+mmscore;
                                    matchindex = {CELL_GAPINY};
                                }}
                            }}else{{
                                if(dp_matrix[p3lt+{CELL_GAPINX}]+mmscore
                                >= dp_matrix[p3lt+{CELL_GAPINY}]+mmscore){{    
                                    matchscore = dp_matrix[p3lt+{CELL_GAPINX}]+mmscore;
                                    matchindex = {CELL_GAPINX};
                                }}else{{
                                    matchscore = dp_matrix[p3lt+{CELL_GAPINY}]+mmscore;
                                    matchindex = {CELL_GAPINY};
                                }}
                            }}

                            if(alignment_type != {ALIGN_LOCAL}){{
                                //emboss の water と needle で優先順位が違うようだ・・・
                                float popenal = (num_cols-col_id-1 == 0)?(start_openal):(openal);
                                float pepenal = (num_cols-col_id-1 == 0)?(start_epenal):(epenal);
                                
                                float qopenal = (num_rows-currentrow-1 == 0)?(start_openal):(openal);
                                float qepenal = (num_rows-currentrow-1 == 0)?(start_epenal):(epenal);
                                
                                if(dp_matrix[p3t+{CELL_MATCH}]+popenal
                                > dp_matrix[p3t+{CELL_GAPINX}]+pepenal){{
                                    if(dp_matrix[p3t+{CELL_MATCH}]+popenal >= dp_matrix[p3t+{CELL_GAPINY}]+popenal){{
                                        gapxscore = dp_matrix[p3t+{CELL_MATCH}]+popenal;
                                        gapxindex = {CELL_MATCH};
                                    }}else{{
                                        gapxscore = dp_matrix[p3t+{CELL_GAPINY}]+popenal;
                                        gapxindex = {CELL_GAPINY};
                                    }}
                                }}else{{
                                    if(dp_matrix[p3t+{CELL_GAPINX}]+pepenal >= dp_matrix[p3t+{CELL_GAPINY}]+popenal){{
                                        gapxscore = dp_matrix[p3t+{CELL_GAPINX}]+pepenal;
                                        gapxindex = {CELL_GAPINX};
                                    }}else{{
                                        gapxscore = dp_matrix[p3t+{CELL_GAPINY}]+popenal;
                                        gapxindex = {CELL_GAPINY};
                                    }}
                                }}
                                if(dp_matrix[p3l+{CELL_MATCH}]+qopenal
                                > dp_matrix[p3l+{CELL_GAPINY}]+qepenal){{
                                    if(dp_matrix[p3l+{CELL_MATCH}]+qopenal >= dp_matrix[p3l+{CELL_GAPINX}]+qopenal){{
                                        gapyscore = dp_matrix[p3l+{CELL_MATCH}]+qopenal;
                                        gapyindex = {CELL_MATCH};
                                    }}else{{
                                        gapyscore = dp_matrix[p3l+{CELL_GAPINX}]+qopenal;
                                        gapyindex = {CELL_GAPINX};
                                    }}
                                }}else{{
                                    if(dp_matrix[p3l+{CELL_GAPINY}]+qepenal >= dp_matrix[p3l+{CELL_GAPINX}]+qopenal){{
                                        gapyscore = dp_matrix[p3l+{CELL_GAPINY}]+qepenal;
                                        gapyindex = {CELL_GAPINY};
                                    }}else{{
                                        gapyscore = dp_matrix[p3l+{CELL_GAPINX}]+qopenal;
                                        gapyindex = {CELL_GAPINX};
                                    }}
                                }}
                            }}else{{
                                if(dp_matrix[p3t+{CELL_MATCH}]+openal
                                >= dp_matrix[p3t+{CELL_GAPINX}]+epenal){{
                                    if(dp_matrix[p3t+{CELL_MATCH}]+openal > dp_matrix[p3t+{CELL_GAPINY}]+openal){{
                                        gapxscore = dp_matrix[p3t+{CELL_MATCH}]+openal;
                                        gapxindex = {CELL_MATCH};
                                    }}else{{
                                        gapxscore = dp_matrix[p3t+{CELL_GAPINY}]+openal;
                                        gapxindex = {CELL_GAPINY};
                                    }}
                                }}else{{
                                    if( dp_matrix[p3t+{CELL_GAPINX}]+epenal > dp_matrix[p3t+{CELL_GAPINY}]+openal){{
                                        gapxscore = dp_matrix[p3t+{CELL_GAPINX}]+epenal;
                                        gapxindex = {CELL_GAPINX};
                                    }}else{{
                                        gapxscore = dp_matrix[p3t+{CELL_GAPINY}]+openal;
                                        gapxindex = {CELL_GAPINY};
                                    }}
                                }}
            
                                if(dp_matrix[p3l+{CELL_MATCH}]+openal
                                >= dp_matrix[p3l+{CELL_GAPINY}]+epenal){{
                                    if(dp_matrix[p3l+{CELL_MATCH}]+openal > dp_matrix[p3l+{CELL_GAPINX}]+openal){{
                                        gapyscore = dp_matrix[p3l+{CELL_MATCH}]+openal;
                                        gapyindex = {CELL_MATCH};
                                    }}else{{
                                        gapyscore = dp_matrix[p3l+{CELL_GAPINX}]+openal;
                                        gapyindex = {CELL_GAPINX};
                                        
                                    }}
            
                                }}else{{
                                    if(dp_matrix[p3l+{CELL_GAPINY}]+epenal > dp_matrix[p3l+{CELL_GAPINX}]+openal){{
                                        gapyscore = dp_matrix[p3l+{CELL_GAPINY}]+epenal;
                                        gapyindex = {CELL_GAPINY};
                                    }}else{{
                                        gapyscore = dp_matrix[p3l+{CELL_GAPINX}]+openal;
                                        gapyindex = {CELL_GAPINX};
                                    }}
                                }}
                                if(matchscore < 0.0){{
                                    matchscore = 0.0;
                                }}
                                if(gapxscore < 0.0){{
                                    gapxscore = 0.0;
                                }}
                                if(gapyscore < 0.0){{
                                    gapyscore = 0.0;
                                }}
                            }}


                            dp_matrix[cp3] = matchscore;
                            dp_matrix[cp3+{CELL_GAPINX}] = gapxscore;
                            dp_matrix[cp3+{CELL_GAPINY}] = gapyscore;

                            if(matchscore > maxscore){{
                                maxpos = currentrow;
                                maxscore = matchscore;
                            }}

                            flag_matrix[currentpos] = (matchindex << 1) + (gapxindex << 3) + (gapyindex << 5) +1;
                            currentrow += 1;
                            
                            prevpos_t = pos_2d_to_1d_rc(currentrow-1,col_id,num_rows,num_cols);
                            prevpos_l = pos_2d_to_1d_rc(currentrow,col_id-1,num_rows,num_cols);
                            prevpos_lt = pos_2d_to_1d_rc(currentrow-1,col_id-1,num_rows,num_cols);
                            currentpos = pos_2d_to_1d_rc(currentrow,col_id,num_rows,num_cols);
                            
                        }}
                    }}
                }}
                max_position[col_id] = maxpos;
                score_buff[col_id] = maxscore;
            }}else{{
                //なんか止まらない？？？
                //barrier(CLK_GLOBAL_MEM_FENCE);
                //if(get_global_id(0) == num_cols-1){{
                if(get_global_id(0) == 0){{
                        backtrack(aligned_seq_a,aligned_seq_b,max_position,score_buff,seqsize,dp_matrix,flag_matrix,alignment_type);
                }}
            }}
        }}
    
        "#,CELL_MATCH=CELL_MATCH
        ,CELL_GAPINX=CELL_GAPINX
        ,CELL_GAPINY=CELL_GAPINY
        ,ALIGN_LOCAL=ALIGN_LOCAL
        ,ALIGN_GLOCAL=ALIGN_GLOCAL
        ,ALIGN_GLOBAL=ALIGN_GLOBAL
        ,POS_TERMINAL=POS_TERMINAL
        );
        let platform = Platform::default();
        let device = Device::first(platform).unwrap_or_else(|e|panic!("{:?}",e));
        let context = Context::builder()
            .platform(platform)
            .devices(device.clone())
            .build().unwrap_or_else(|e|panic!("{:?}",e));
        let program = Program::builder()
            .devices(device)
            .src(src_test)
            .build(&context).unwrap_or_else(|e|panic!("{:?}",e));

        let queue = Queue::new(&context, device, None).unwrap_or_else(|e|panic!("{:?}",e));
        let matrix_size:usize = (max_length+1)*(max_length+1);
        let dp_matrix:VecAndBuffer<f32> = VecAndBuffer::new(vec![0.0_f32;matrix_size*3],&queue);
        let flag_matrix:VecAndBuffer<u8> = VecAndBuffer::new(vec![0_u8;matrix_size],&queue);
        let seq_length:VecAndBuffer<i32> = VecAndBuffer::new(vec![0_i32;2],&queue);
        let kernel_control:VecAndBuffer<u8> = VecAndBuffer::new(vec![0_u8;1],&queue);

        let mut pgo = go;
        let mut pge = ge;
        if pgo > 0.0{
            pgo *= -1.0_f32;
        }
        if pge > 0.0{
            pge *= -1.0_f32;
        }
        let mut spgo = 0.0;
        let mut spge = 0.0;
        if alignment_type == ALIGN_GLOBAL{
            spgo = pgo;
            spge = pge;
        }

        let scoring_matrix_vec_ = ss.get_vec_score();

        let scoring_matrix_vec =  match 
        ocl::Buffer::<f32>::builder()
        .queue(queue.clone())
        .flags(ocl::MemFlags::new().read_write().copy_host_ptr())
        .len(scoring_matrix_vec_.len())
        .copy_host_slice(scoring_matrix_vec_.as_slice())
        .build(){
            Ok(x) => {x},
            Err(e) =>{
                panic!("{:?}",e);
            }
        };
        
        scoring_matrix_vec.write(&scoring_matrix_vec_).enq().unwrap_or_else(|e|{panic!("{:?}",e)});
        let seq_a:VecAndBuffer<i32> = VecAndBuffer::new(vec![0_i32;max_length],&queue);
        let seq_b:VecAndBuffer<i32> = VecAndBuffer::new(vec![0_i32;max_length],&queue);
        
        let aligned_seq_a:VecAndBuffer<i32> = VecAndBuffer::new(vec![POS_TERMINAL;max_length*2+1],&queue);
        let max_position:VecAndBuffer<i32> = VecAndBuffer::new(vec![POS_TERMINAL;max_length],&queue);
        let aligned_seq_b:VecAndBuffer<i32> = VecAndBuffer::new(vec![POS_TERMINAL;max_length*2+1],&queue);
        let score_buff:VecAndBuffer<f32> = VecAndBuffer::new(vec![0.0;max_length],&queue);

        let kernel_fill = Kernel::builder()
        .program(&program)
        .name("fill_matrix")
        .queue(queue.clone())
        .global_work_size(seq_a.vector.len()+1)
        .arg(&seq_a.buffer)
        .arg(&seq_b.buffer)
        .arg(&seq_length.buffer)
        .arg(&scoring_matrix_vec)
        .arg(ss.get_num_columns() as i32)
        .arg(&dp_matrix.buffer)
        .arg(&flag_matrix.buffer)
        .arg(spgo as f32)
        .arg(spge as f32)
        .arg(pgo as f32)
        .arg(pge as f32)
        .arg(&aligned_seq_a.buffer)
        .arg(&aligned_seq_b.buffer)
        .arg(&max_position.buffer)
        .arg(&score_buff.buffer)
        .arg(alignment_type as i32)
        .arg(&kernel_control.buffer)
        .build().unwrap_or_else(|e|panic!("{:?}",e));

        return OpenCLSequenceAlignment{
            scoring_matrix:ss
            ,scoring_matrix_vec:scoring_matrix_vec
            ,o_penalty:pgo
            ,e_penalty:pge
            ,kernel_fill:kernel_fill
            ,start_openal:spgo
            ,start_epenal:spge
            ,seq_a:seq_a
            ,seq_b:seq_b
            ,aligned_seq_a:aligned_seq_a
            ,aligned_seq_b:aligned_seq_b
            ,max_position:max_position
            ,score_buff:score_buff
            ,seq_length:seq_length
            ,num_cols:0
            ,num_rows:0
            ,alignment_type:alignment_type
            ,ocl_platform:platform
            ,ocl_device:device
            ,ocl_context:context
            ,ocl_program:program
            ,ocl_queue:queue
            ,dp_matrix:dp_matrix
            ,flag_matrix:flag_matrix
            ,kernel_control:kernel_control
        };
    }

    pub fn clear_buffer(&mut self){
        //self.cells = vec![vec![]];
    }

    pub fn set_scoring_matrix(&mut self,scoringmatrix:Box<dyn ScoringMatrix>){
        self.scoring_matrix = scoringmatrix;
    }
    pub fn prepare(&mut self,s1:&SeqData,s2:&SeqData){
        
        let buffremake_diff:i64 = 1000;

        let seq_a = self.scoring_matrix.seq_to_index(&s1,None);
        let seq_b = self.scoring_matrix.seq_to_index(&s2,None);
        self.num_cols = seq_a.len()+1;
        self.num_rows = seq_b.len()+1;
        self.seq_length.vector[0] = seq_a.len() as i32;
        self.seq_length.vector[1] = seq_b.len() as i32;
        self.seq_length.write();
        
        if self.num_cols*self.num_rows*3 > 2147483647{
            panic!("I think the sequences are too long!");
        }

        let mut buffer_updated:bool = false;
        if (seq_a.len() > self.seq_a.vector.len()) || (seq_a.len() as i64) < (self.seq_a.vector.len() as i64) - buffremake_diff {
            self.seq_a = VecAndBuffer::new(vec![0_i32;seq_a.len()+((buffremake_diff/4) as usize)],&self.ocl_queue);
            self.max_position = VecAndBuffer::new(vec![0_i32;seq_a.len()+((buffremake_diff/4) as usize)],&self.ocl_queue);
            self.score_buff = VecAndBuffer::new(vec![0_f32;seq_a.len()+((buffremake_diff/4) as usize)],&self.ocl_queue);
          
            buffer_updated = true;
        }
        
        if seq_b.len() > self.seq_b.vector.len() || (seq_b.len() as i64) < (self.seq_b.vector.len() as i64) - buffremake_diff{
            self.seq_b = VecAndBuffer::new(vec![0_i32;seq_b.len()+((buffremake_diff/4) as usize)],&self.ocl_queue);
            buffer_updated = true;
        }
        if buffer_updated{
            let llen = self.seq_a.vector.len()+self.seq_b.vector.len()+2;
            self.aligned_seq_a = VecAndBuffer::new(vec![POS_TERMINAL;llen],&self.ocl_queue);
            self.aligned_seq_b = VecAndBuffer::new(vec![POS_TERMINAL;llen],&self.ocl_queue);
        }
        for ii in 0..seq_a.len(){
            self.seq_a.vector[ii] = seq_a[ii];
        }
        for ii in 0..seq_b.len(){
            self.seq_b.vector[ii] = seq_b[ii];
        }
        
        self.seq_a.write();
        self.seq_b.write();
        self.aligned_seq_a.write();
        self.aligned_seq_b.write();
        self.score_buff.vector[0] = 0.0;
        self.score_buff.write();


        let matrix_size:usize = self.num_cols*self.num_rows;
        if matrix_size > self.flag_matrix.vector.len() || (matrix_size as i64) < (self.flag_matrix.vector.len() as i64) - buffremake_diff*buffremake_diff{
            self.dp_matrix = VecAndBuffer::new(vec![0.0_f32;(matrix_size+((buffremake_diff*buffremake_diff/16) as usize))*3],&self.ocl_queue);
            self.flag_matrix = VecAndBuffer::new(vec![0_u8;matrix_size+((buffremake_diff*buffremake_diff/16) as usize)],&self.ocl_queue);
            buffer_updated = true;
        }


        if buffer_updated{
            self.kernel_fill = Kernel::builder()
            .program(&self.ocl_program)
            .name("fill_matrix")
            .queue(self.ocl_queue.clone())
            .global_work_size(self.num_cols)
            .arg(&self.seq_a.buffer)
            .arg(&self.seq_b.buffer)
            .arg(&self.seq_length.buffer)
            .arg(&self.scoring_matrix_vec)
            .arg(self.scoring_matrix.get_num_columns() as i32)
            .arg(&self.dp_matrix.buffer)
            .arg(&self.flag_matrix.buffer)
            .arg(self.start_openal as f32)
            .arg(self.start_epenal as f32)
            .arg(self.o_penalty as f32)
            .arg(self.e_penalty as f32)
            .arg(&self.aligned_seq_a.buffer)
            .arg(&self.aligned_seq_b.buffer)
            .arg(&self.max_position.buffer)
            .arg(&self.score_buff.buffer)
            .arg(self.alignment_type as i32)
            .arg(&self.kernel_control.buffer)
            .build().unwrap_or_else(|e|panic!("{:?}",e));
        }
    }
    
    pub fn pos_2d_to_1d_rc_(&self,r:usize,c:usize)->usize{
        return r*self.num_cols+c;
    }
    

    pub fn pos_2d_to_1d_rc(&self,r_:usize,c_:usize)->usize{
        let r = r_ as i64;
        let c = c_ as i64;
        let k = (c +r -self.num_rows as i64+1).max(0);
        let l = (c +r -self.num_cols as i64+1).max(0);
        let ks = k*self.num_rows as i64;
        let ls = l*self.num_cols as i64;
        let pt = l*k;
        let ps = ks+ls-pt;
        let cc = c-k;
        let rr = r-l;
        let zt = (cc+rr)*(cc+rr+1)/2+rr;

        return (ps+zt) as usize;
    }

    pub fn backtrack(&mut self)->(Vec<i64>,Vec<i64>,f32){
        let mut ret1:Vec<i64> = vec![];
        let mut ret2:Vec<i64> = vec![];
        self.aligned_seq_a.read();
        self.aligned_seq_b.read();
        self.score_buff.read();
        let maxscore = self.score_buff.vector[0];
        
        let mut ppos:usize = 0;
        loop{
            if self.aligned_seq_a.vector[ppos] != POS_TERMINAL{
                ret1.push(self.aligned_seq_a.vector[ppos] as i64);
            }else{
                break;
            }
            ppos += 1;
        }//同じ長さになるはずではある
        
        let mut ppos:usize = 0;
        loop{
            if self.aligned_seq_b.vector[ppos] != POS_TERMINAL{
                ret2.push(self.aligned_seq_b.vector[ppos] as i64);
            }else{
                break;
            }
            ppos += 1;
        }
        ret1.reverse();
        ret2.reverse();
        return (ret1,ret2,maxscore);
    }
    pub fn align(&mut self,s1:&SeqData,s2:&SeqData,retain_all:bool)->(Vec<String>,Vec<String>,f32){
        self.prepare(&s1, &s2);
        
        unsafe {
            self.kernel_control.vector[0] = 0;
            self.kernel_control.write();
            self.kernel_fill.cmd()
            .queue(&self.ocl_queue)
            .global_work_offset(ocl::SpatialDims::new(Some(0),None, None).unwrap())
            .global_work_size(self.num_cols)
            .local_work_size(self.kernel_fill.default_local_work_size())
            .enq().unwrap_or_else(|e|panic!("{:?}",e));
            
            self.kernel_control.vector[0] = 1;
            self.kernel_control.write();
            self.kernel_fill.cmd()
            .queue(&self.ocl_queue)
            .global_work_offset(ocl::SpatialDims::new(Some(0),None, None).unwrap())
            .global_work_size(self.num_cols)
            .local_work_size(self.kernel_fill.default_local_work_size())
            .enq().unwrap_or_else(|e|panic!("{:?}",e));
            
            self.kernel_control.vector[0] = 2;
            self.kernel_control.write();
            self.kernel_fill.cmd()
            .queue(&self.ocl_queue)
            .global_work_offset(ocl::SpatialDims::new(Some(0),None, None).unwrap())
            .global_work_size(1)
            .local_work_size(self.kernel_fill.default_local_work_size())
            .enq().unwrap_or_else(|e|panic!("{:?}",e));
            
        }
        
        //self.flag_matrix.read();
        //self.dp_matrix.read();


        //print!("+++++++++++++++++++++++++++++++++++++++++++++++++{:?}",self.flag_matrix.vector);
        let res = self.backtrack();
        //print!("+++++{:?}",res);
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



pub trait ScoringMatrix{
    fn get_score(&self,a:usize,b:usize)->f32;
    fn get_score_str(&self,a:&str,b:&str)->f32;
    fn seq_to_index(&self,ss:&SeqData,partial_region:Option<usize>)->Vec<i32>;
    fn set_score(&mut self,a:usize,b:usize,s:f32);
    fn prepare(&mut self,a:&SeqData,b:&SeqData);
    fn get_vec_score(&self)->Vec<f32>;
    fn get_num_columns(&self)->usize;
}
/* not tested
#[allow(dead_code)]
pub struct PositionSpecificMatrix{
    pub scores:Vec<f32>,
    pub a_length:usize,
    pub b_length:usize
}

impl ScoringMatrix for PositionSpecificMatrix{
    fn get_vec_score(&self)->Vec<f32>{
        return self.scores.clone();
    }
    fn get_num_columns(&self) ->usize {
        return self.b_length;
    }
    fn get_score_str(&self,_a:&str,_b:&str)->f32{
        panic!("not implemented");
    }
    fn get_score(&self,a:usize,b:usize)->f32{
        return self.scores[a+b*self.a_length];
    }
    fn seq_to_index(&self,ss:&SeqData,partial_region:Option<usize>)->Vec<i32>{
        if let Some(x) = partial_region{
        return (0..x).into_iter().map(|m| m as i32).collect();
        }else{
            return (0..ss.seq.len()).into_iter().map(|m| m as i32).collect();
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
*/
#[allow(dead_code)]
pub struct SubstitutionMatrix{
    pub string_to_index:HashMap<String,usize> ,
    pub index_to_string:Vec<String>,
    pub scores:Vec<Vec<f32>>,
    
}
impl ScoringMatrix for SubstitutionMatrix{
    fn get_vec_score(&self)->Vec<f32>{
        let slen:usize = self.scores.len();
        let mut ret:Vec<f32> = vec![0.0;slen*slen];
        for rr in 0..slen{
            for cc in 0..slen{
                ret[rr*slen+cc] = self.scores[rr][cc];
            } 
        }
        return ret;
    }
    fn get_num_columns(&self)->usize{
        return self.scores[0].len();
    }


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
    fn seq_to_index(&self,ss:&SeqData,partial_region:Option<usize>)->Vec<i32>{
        if let Some(x) = partial_region{
            let mut ret:Vec<i32> = vec![];
            for xx in 0..x{
                ret.push(self.get_string_index(&ss.seq[xx]) as i32);
            }
            return ret;
        }else{ 
            return ss.seq.iter().map(|m| self.get_string_index(m) as i32).collect();
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
    
    //let seq1_ = "AAAAAAAAA".to_string();
    //let seq2_ = "AAAAAAAAA".to_string();

    let mut sw = OpenCLSequenceAlignment::new(100,Box::new(sm),10.0,0.5,ALIGN_LOCAL);

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
    let mut sw = OpenCLSequenceAlignment::new(100,Box::new(sm),10.0,0.5,ALIGN_GLOCAL);
    let res = sw.align(&seq1,&seq2,true);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    
    assert_eq!(res.2,51.5);
    assert_eq!(r1,"----CATTAGATGACT-----GAAAGCAAG----------TACTGGTC------TCTTAAACCATTTAATAGTAAATTAGCACTTACTTCTAATGA");
    assert_eq!(r2,"ACTTCTCTAGCTCAGTTGGTAGAGCGCAAGGCTTTTAACCTTGTGGTCGTGGGTTC--AAACCCCATGATGG-------GCA--------------");

    let sm = SubstitutionMatrix::get_mat_matrix(5.0_f32,-4.0_f32);
    let mut sw = OpenCLSequenceAlignment::new(100,Box::new(sm),10.0,0.5,ALIGN_GLOBAL);
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
    let mut sw = OpenCLSequenceAlignment::new(100,Box::new(sm),10.0,0.5,ALIGN_GLOCAL);
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
    let mut sw = OpenCLSequenceAlignment::new(100,Box::new(sm),10.0,0.5,ALIGN_GLOBAL);
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

    
    let mut sm = SubstitutionMatrix::get_blosum62_matrix();
    for vv in sm.scores.iter_mut(){
        for vvv in vv.iter_mut(){
            *vvv *= 10.0; 
        }
    }
    
    let mut sw = OpenCLSequenceAlignment::new(100,Box::new(sm),100.0,5.0,ALIGN_LOCAL);
    let seq1_ = "MAGELTPEEEAQYKKAFSAVDTDGNGTINAQELGAALKATGKNLSEAQLRKLISEVDSDGDGEISFQEFLTAAKKARAGLEDLQVAFRAFDQDGDGHITVDELRRAMAGLGQPLPQEELDAMIREADVDQDGRVNYEEFARMLAQE".to_owned();
    let seq2_ = "MNFTPTHTPVCRKRTVVSKRGVAVSGPTKRRGMADSLESTPLPSPEDRLAKLHPSKELLEYYQKKMAECEAENEDLLKKLELYKEACEGQHKLECDLQQREEEIAELQKALSDMQVCLFQEREHVLRLYSENDRLRIRELEDKKKIQNLLALVGTDAGEVTYFCKEPPHKVTILQKTIQAVGECEQSESSAFKADPKISKRRPSRERKESSEHYQRDIQTLILQVEALQAQLGEQTKLSREQIEGLIEDRRIHLEEIQVQHQRNQNKIKELTKNLHHTQELLYESTKDFLQLRSENQNKEKSWMLEKDNLMSKIKQYRVQCKKKEDKIGKVLPVMHESHHAQSEYIKSLKDKLVQEKKLSNMYQEQCISLEEELARIREEEGMRREIFKDRTNKMGKRLQIMTKRYEALERRRILEVEGFKTDIKVLRQKLKDLEQMLYKATVNARANQDLALLCEVRDSNRRAHKIQGELKNLKSKVFGLENELRLC".to_owned();

    let seq1 = SeqData::create("".to_string(),"".to_string(),seq1_);
    let seq2 = SeqData::create("".to_string(),"".to_string(),seq2_);
    //let sm = ScoringMatrix::get_blosum62_matrix();
    //let mut sw = SequenceAlignment::new(sm,10.0,0.5,ALIGN_GLOBAL);
    let res = sw.align(&seq1,&seq2,false);
    let r1 = res.0.iter().fold("".to_string(),|s,m|s+m);
    let r2 = res.1.iter().fold("".to_string(),|s,m|s+m);
    assert_eq!(r1,"EEEAQYKKAFSAVDTDGNGTINAQELGAALKATGKNLSEAQLRKLISEVDSDGDGEISFQEFLTAAKKARAGLEDL-QVAFRA");
    assert_eq!(r2,"EEEGMRREIFKD---------RTNKMGKRLQIMTKRY-EALERRRILEVEG----------FKTDIKVLRQKLKDLEQMLYKA");

}
