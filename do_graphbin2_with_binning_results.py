#Working on the folder contain '04.assembly' and '06.binning' folder generated from metapi workflow: https://github.com/ohmeta/metapi
import os
import time
import datetime
import sys

sample_num=[0,2,3,5,6,21,22,24,25]
binners=['concoct','maxbin2','metabat2','vamb']
for sample_num in sample_num:
    binning_out_dir=f'CAMI_UT_{sample_num}.metaspades.out'
    text=f'**********\nStarting graphbin2 with {binning_out_dir}/ ...\n'
    for words in text:
        sys.stdout.write(words)
        sys.stdout.flush()
        time.sleep(0.01)
    
    os.system(f'cp 04.assembly/scaftigs/CAMI_UT_{sample_num}.metaspades.out/CAMI_UT_{sample_num}.metaspades.scaftigs.* .;\
    gunzip *.gz; mv CAMI_UT_{sample_num}.metaspades.scaftigs.fa CAMI_UT_{sample_num}.metaspades.scaftigs.fasta')
    
    with open(f'CAMI_UT_{sample_num}.metaspades.scaftigs.fasta',"r") as f:
        fasta_file=f.readlines()
    contigs={}
    i=0
    while i<len(fasta_file):
        if fasta_file[i][0]=='>':
            contig_name=fasta_file[i].rstrip().lstrip('>')
            contigs[contig_name]=[]
            contigs[contig_name].append(fasta_file[i+1])
            i=i+2
        else:
            contigs[contig_name].append(fasta_file[i])
            i=i+1
                
    for binner in binners: #Mustn't be 'for binner in binner' cause its the second loop
        graphbin_out_dir=f'{binning_out_dir}/{binner}_graphbin2'
        text=f'Processing graphbin2 of {graphbin_out_dir} ...\n'
        for words in text:
            sys.stdout.write(words)
            sys.stdout.flush()
            time.sleep(0.01)
        start=datetime.datetime.now()
        
        os.system(f'mkdir -p {graphbin_out_dir}')
        os.system(f'python ../GraphBin2/support/prepResult.py --binned 06.binning/bins/CAMI_UT_{sample_num}.metaspades.out/{binner}/ --output .')
        if binner=='vamb':# Cut the "Snc" words in the vamb generated mag's contig name, '>SncNode_123_length_123_coverage_1.23' => '>Node_123_length_123_coverage_1.23'
            os.system("s=$(sed -n '1p' initial_contig_bins.csv | awk -F 'NODE' '{print $1}' | cut -d '>' -f2);sed -i \"s/$s//g\" initial_contig_bins.csv")
        os.system(f'python ../GraphBin2/graphbin2 --assembler spades --contigs CAMI_UT_{sample_num}.metaspades.scaftigs.fasta \
        --graph CAMI_UT_{sample_num}.metaspades.scaftigs.gfa --paths CAMI_UT_{sample_num}.metaspades.scaftigs.paths \
        --binned initial_contig_bins.csv --nthreads 20 --output .')
        
        bin_ids=[]
        with open("bin_ids.csv","r") as f:
            bin_ids_file=f.readlines()
        for i in range(len(bin_ids_file)):
            bin_ids.append(bin_ids_file[i].split(',')[0])

        graphbin_fa_file=[]
        for i in bin_ids: #CAMI_UT_0.metaspades.maxbin2.bin.5.fa
            j=i.split('.') #['CAMI_UT_0', 'metaspades', 'maxbin2', 'bin', '5', 'fa']
            a='.'.join(j[0:3])
            b='.'.join(j[3:])
            c=f'{a}.graphbin2.{b}' #CAMI_UT_0.metaspades.maxbin2.graphbin2.bin.5.fa
            graphbin_fa_file.append(c)
            os.system('touch %s' % c) 
            
        graphbin_out_csv=[]
        with open("graphbin2_output.csv","r") as f:
            graphbin_out_csv=f.readlines() #'NODE_1_length_2618738_cov_87.904107,2\n'
        for i in graphbin_out_csv:
            sentence=i.rstrip().split(',')
            contig_name=sentence[0]
            fa_file_num=int(sentence[1])-1
            with open(graphbin_fa_file[fa_file_num],"a") as f:
                f.writelines('>'+contig_name+'\n')
                f.writelines(contigs[contig_name])
                
        f=os.popen("rg Found graphbin2.log | cut -d ' ' -f6-") #Found 20 multi-labelled contigs
        print(f.read()+f'Check {graphbin_out_dir}/graphbin2.log for more details.') #Check {binning_out_dir}/{binner}_graphbin2/graphbin2.log for more details.
        
        os.system(f'mv *.fa {graphbin_out_dir}; mv *.csv {graphbin_out_dir}; mv *.log {graphbin_out_dir}')
        end=datetime.datetime.now()
        running_time=end-start
        print(f'{graphbin_out_dir} finished. Running time: {running_time}\n')
        
    dastool_prefix=f'CAMI_UT_{sample_num}.metaspades.dastools.bin'
    os.system(f'mkdir -p {binning_out_dir}/dastool')
    text=f'\nStarting DAS_Tool with {binning_out_dir}/ ...\n'
    for words in text:
        sys.stdout.write(words)
        sys.stdout.flush()
        time.sleep(0.01)
    start=datetime.datetime.now()

    for binner in binners:
        graphbin_out_dir=f'{binning_out_dir}/{binner}_graphbin2'
        os.system(f'Fasta_to_Scaffolds2Bin.sh -i {graphbin_out_dir} -e fa > {dastool_prefix}.{binner}_graphbin2_scaffolds2bin.csv')

    with open(f'{dastool_prefix}.sh','w') as f:
        f.write(f'#!/bin/bash\n\
                DAS_Tool \\\n\
                --bins {dastool_prefix}.{binners[0]}_graphbin2_scaffolds2bin.csv,{dastool_prefix}.{binners[1]}_graphbin2_scaffolds2bin.csv,{dastool_prefix}.{binners[2]}_graphbin2_scaffolds2bin.csv,{dastool_prefix}.{binners[3]}_graphbin2_scaffolds2bin.csv \\\n\
                --contigs CAMI_UT_{sample_num}.metaspades.scaftigs.fasta \\\n\
                --labels {binners[0]},{binners[1]},{binners[2],{binners[3]} \\\n\
                --outputbasename {dastool_prefix} \\\n\
                --search_engine diamond \\\n\
                --write_bins 1 \\\n\
                --threads 20')
    os.system("chmod +x *.sh; bash *.sh")
    os.system(f'mv {dastool_prefix}* {binning_out_dir}/dastool')

    end=datetime.datetime.now()
    running_time=end-start
    print(f'DAS_Tool with {binning_out_dir}/ finished. Running time: {running_time}\nChecking {binning_out_dir}/dastool for more details.\n')
        
    os.system(f'rm CAMI_UT_{sample_num}.metaspades.scaftigs.*')
    print(f'{binning_out_dir}/ finished.\n**********\n')
print("All processes finished!")
