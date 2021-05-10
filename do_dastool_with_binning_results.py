#Working on the folder with '04.assembly' and '06.binning' folder generated from metapi workflow: https://github.com/ohmeta/metapi, and you also need to generate the relative sample folder like CAMI_UT_{sample_num}.metaspades.out beforehand!

##initialize directory structure
#bash -c 'for i in $(ls ../concoct_maxbin2_metabat2_vamb_graphbin2_dastool/ | cut -d "_" -f3 | cut -d "." -f1);do mkdir CAMI_UT_$i.metaspades.out; cp -r 06.binning/bins/CAMI_UT_$i.metaspades.out/concoct/ CAMI_UT_$i.metaspades.out/;cp -r 06.binning/bins/CAMI_UT_$i.metaspades.out/maxbin2/ CAMI_UT_$i.metaspades.out/;cp -r 06.binning/bins/CAMI_UT_$i.metaspades.out/metabat2/ CAMI_UT_$i.metaspades.out/;cp -r 06.binning/bins/CAMI_UT_$i.metaspades.out/vamb/ CAMI_UT_$i.metaspades.out/;done'

import os
import time
import datetime
import sys

sample_num=[0,2,3,5,6,21,22,24,25]
binners=['concoct','maxbin2','metabat2','vamb']
for sample_num in sample_num:
    binning_out_dir=f'CAMI_UT_{sample_num}.metaspades.out'
    dastool_prefix=f'CAMI_UT_{sample_num}.metaspades.dastools.bin'
    os.system(f'mkdir -p {binning_out_dir}/dastool')
    text=f'\nStarting DAS_Tool with {binning_out_dir}/ ...\n'
    for words in text:
        sys.stdout.write(words)
        sys.stdout.flush()
        time.sleep(0.01)
    start=datetime.datetime.now()
    
    os.system(f'cp 04.assembly/scaftigs/CAMI_UT_{sample_num}.metaspades.out/CAMI_UT_{sample_num}.metaspades.scaftigs.fa.gz .;\
    gunzip *.gz; mv CAMI_UT_{sample_num}.metaspades.scaftigs.fa CAMI_UT_{sample_num}.metaspades.scaftigs.fasta')
    
    for binner in binners:
        graphbin_out_dir=f'{binning_out_dir}/{binner}_graphbin2'
        os.system(f'Fasta_to_Scaffolds2Bin.sh -i {graphbin_out_dir} -e fa > {dastool_prefix}.{binner}_graphbin2_scaffolds2bin.csv')
    
    with open(f'{dastool_prefix}.sh','w') as f:
        f.write(f'#!/bin/bash\n\
                DAS_Tool \\\n\
                --bins {dastool_prefix}.{binners[0]}_graphbin2_scaffolds2bin.csv,{dastool_prefix}.{binners[1]}_graphbin2_scaffolds2bin.csv,{dastool_prefix}.{binners[2]}_graphbin2_scaffolds2bin.csv \\\n\
                --contigs CAMI_UT_{sample_num}.metaspades.scaftigs.fasta \\\n\
                --labels {binners[0]},{binners[1]},{binners[2]} \\\n\
                --outputbasename {dastool_prefix} \\\n\
                --search_engine diamond \\\n\
                --write_bins 1 \\\n\
                --threads 20')
    os.system("chmod +x *.sh; bash *.sh")
    os.system(f'mv {dastool_prefix}* {binning_out_dir}/dastool')
    
    end=datetime.datetime.now()
    running_time=end-start
    print(f'DAS_Tool with {binning_out_dir}/ finished. Running time: {running_time}\nChecking {binning_out_dir}/dastool for more details.\n')
print("All processes finished!")
