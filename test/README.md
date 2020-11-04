1-4 can be skipped.

==== Please use 1-4 when you want to construct more than 100 example files. 100 is not enough to cover all cases. I'm too stupid to make a perfect testset....====
1. Install emboss package, perl, python(3), java, and cargo(rust).
2. Prepare a multiple fasta file. I used swiss-prot entries of Homo sapiens.
3. Change $sourcefas variable in embpss_run.pl to the file path you made in step 2. & $embossneedle, $embosswater, $num_trials, can manipulate other settings.
4. Type `perl emboss_run.pl`.

 ** As it seems that EMBOSS needle has a bug, test may fail with files produced with this process. I reported but got server unreachable error.

====

5. Type `perl run_java.pl`, `perl run_py.pl`, `perl run_pl.pl`, `perl run_rust.pl`.
6. Type `perl check_results.pl`.
7. If default dataset was used, 

OK
Checked 600 results.

will be shown as message.

\* Only rust version produces global and glocal (global + endweight==0) alignments.
