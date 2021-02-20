# smithwaterman

public domain smith-waterman (& some of them have needleman-wunsch) algorithm scripts which can be used as a module in your own applications.  
  
When benchmarked with 1(256aa) vs 5000 (<= 3000aa) swiss-prot(downloaded in 2021/02/21) proteins, 
Rust+OpenCL version is about 2.5 times SLOWER than the Rust (cpu) version...(37 sec vs 15 sec)🙄  
  
When benchmarked with 1(1165aa) vs 5000 (1000aa <= len <= 3000aa) swiss-prot(downloaded in 2021/02/21) proteins, Rust+OpenCL version is about 1.2 times faster than the Rust (cpu) version...(220 sec vs 262 sec) 🙄  

