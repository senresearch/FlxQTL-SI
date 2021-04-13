## implemented in linux machine
#scans :3 traits
(time gemma.3 -g ../testdata/mouse_hs1940_genochr1.txt -p ../testdata/mouse_hs1940_imp_pheno.txt -n 1 2 3 -k ./output/kinship.cXX.txt -lmm 2 -o lrt3_ch1-v3 -no-check)2>> gemma_runtime3trait.txt

 (time gemma.3 -g ../testdata/mouse_hs1940_genochr3.txt -p ../testdata/mouse_hs1940_imp_pheno.txt -n 1 2 3 -k ./output/kinship.cXX.txt -lmm 2 -o lrt3_ch3-v3 -no-check)2>> gemma_runtime3trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940_genochr7.txt -p ../testdata/mouse_hs1940_imp_pheno.txt -n 1 2 3 -k ./output/kinship.cXX.txt -lmm 2 -o lrt3_ch7-v3 -no-check)2>> gemma_runtime3trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940_genochr11.txt -p ../testdata/mouse_hs1940_imp_pheno.txt -n 1 2 3 -k ./output/kinship.cXX.txt -lmm 2 -o lrt3_ch11-v3 -no-check)2>> gemma_runtime3trait.txt

 (time gemma.3 -g ../testdata/mouse_hs1940_genochr15.txt -p ../testdata/mouse_hs1940_imp_pheno.txt -n 1 2 3 -k ./output/kinship.cXX.txt -lmm 2 -o lrt3_ch15-v3 -no-check)2>> gemma_runtime3trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_imp_pheno.txt -n 1 2 3 -k ./output/kinship.cXX.txt -lmm 2 -o lrt3_ch19-v3 -no-check)2>> gemma_runtime3trait.txt

#scans :6 traits
(time gemma.3 -g ../testdata/mouse_hs1940_genochr1.txt -p ../testdata/mouse_hs1940_simulated6traits.txt -n 1 2 3 4 5 6 -k ./output/kinship.cXX.txt -lmm 2 -o lrt6_ch1-v3 -no-check)2>> gemma_runtime6trait.txt

 (time gemma.3 -g ../testdata/mouse_hs1940_genochr3.txt -p ../testdata/mouse_hs1940_simulated6traits.txt -n 1 2 3 4 5 6 -k ./output/kinship.cXX.txt -lmm 2 -o lrt6_ch3-v3 -no-check)2>> gemma_runtime6trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940_genochr7.txt -p ../testdata/mouse_hs1940_simulated6traits.txt -n 1 2 3 4 5 6 -k ./output/kinship.cXX.txt -lmm 2 -o lrt6_ch7-v3 -no-check)2>> gemma_runtime6trait.txt

 (time gemma.3 -g ../testdata/mouse_hs1940_genochr11.txt -p ../testdata/mouse_hs1940_simulated6traits.txt -n 1 2 3 4 5 6 -k ./output/kinship.cXX.txt -lmm 2 -o lrt6_ch11-v3 -no-check)2>> gemma_runtime6trait.txt

 (time gemma.3 -g ../testdata/mouse_hs1940_genochr15.txt -p ../testdata/mouse_hs1940_simulated6traits.txt -n 1 2 3 4 5 6 -k ./output/kinship.cXX.txt -lmm 2 -o lrt6_ch15-v3 -no-check)2>> gemma_runtime6trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_simulated6traits.txt -n 1 2 3 4 5 6 -k ./output/kinship.cXX.txt -lmm 2 -o lrt6_ch19-v3 -no-check)2>> gemma_runtime6trait.txt

#scans :9 traits
(time gemma.3 -g ../testdata/mouse_hs1940_genochr1.txt -p ../testdata/mouse_hs1940_simulated9traits.txt -n 1 2 3 4 5 6 7 8 9 -k ./output/kinship.cXX.txt -lmm 2 -o lrt9_ch1-v3 -no-check)2>> gemma_runtime9trait.txt

 (time gemma.3 -g ../testdata/mouse_hs1940_genochr3.txt -p ../testdata/mouse_hs1940_simulated9traits.txt -n 1 2 3 4 5 6 7 8 9 -k ./output/kinship.cXX.txt -lmm 2 -o lrt9_ch3-v3 -no-check)2>> gemma_runtime9trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940_genochr7.txt -p ../testdata/mouse_hs1940_simulated9traits.txt -n 1 2 3 4 5 6 7 8 9 -k ./output/kinship.cXX.txt -lmm 2 -o lrt9_ch7-v3 -no-check)2>> gemma_runtime9trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940_genochr11.txt -p ../testdata/mouse_hs1940_simulated9traits.txt -n 1 2 3 4 5 6 7 8 9 -k ./output/kinship.cXX.txt -lmm 2 -o lrt9_ch11-v3 -no-check)2>> gemma_runtime9trait.txt

 (time gemma.3 -g ../testdata/mouse_hs1940_genochr15.txt -p ../testdata/mouse_hs1940_simulated9traits.txt -n 1 2 3 4 5 6 7 8 9 -k ./output/kinship.cXX.txt -lmm 2 -o lrt9_ch15-v3 -no-check)2>> gemma_runtime9trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_simulated9traits.txt -n 1 2 3 4 5 6 7 8 9 -k ./output/kinship.cXX.txt -lmm 2 -o lrt9_ch19-v3 -no-check)2>> gemma_runtime9trait.txt

#scans;12 traits
(time gemma.3 -g ../testdata/mouse_hs1940_genochr1.txt -p ../testdata/mouse_hs1940_simulated12traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 -k ./output/kinship.cXX.txt -lmm 2 -o lrt12_ch1-v3 -no-check)2>> gemma_runtime12trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940_genochr3.txt -p ../testdata/mouse_hs1940_simulated12traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 -k ./output/kinship.cXX.txt -lmm 2 -o lrt12_ch3-v3 -no-check)2>> gemma_runtime12trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940_genochr7.txt -p ../testdata/mouse_hs1940_simulated12traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 -k ./output/kinship.cXX.txt -lmm 2 -o lrt12_ch7-v3 -no-check)2>> gemma_runtime12trait.txt


(time gemma.3 -g ../testdata/mouse_hs1940_genochr11.txt -p ../testdata/mouse_hs1940_simulated12traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 -k ./output/kinship.cXX.txt -lmm 2 -o lrt12_ch11-v3 -no-check)2>> gemma_runtime12trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940_genochr15.txt -p ../testdata/mouse_hs1940_simulated12traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 -k ./output/kinship.cXX.txt -lmm 2 -o lrt12_ch15-v3 -no-check)2>> gemma_runtime12trait.txt

(time gemma.3 -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_simulated12traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 -k ./output/kinship.cXX.txt -lmm 2 -o lrt12_ch19-v3 -no-check)2>> gemma_runtime12trait.txt