# smithwaterman

public domain smith-waterman (& some of them have needleman-wunsch) algorithm scripts which can be used as a module in your own applications.  
  
Rust+OpenCL version is about 2.5 times SLOWER than the Rust (cpu) version...🙄  
(Benchmarked with 1(256aa) vs 5000 (<= 3000aa) swiss-prot(downloaded in 2021/02/21) proteins..(picked up from the head). 1 vs 1 x5000 would be much much SLOWER.)    

