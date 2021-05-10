import os
import time
import datetime
import sys
import pprint
import pickle

def chose_candidate_or_perfect_cluster(path_to_perfect_cluster=None,path_to_candidate_cluster=None,\
            path_to_checkm_result=None,path_to_mag=None,show_PC=False,show_CC=False,cluster_type=None,cluster_num=None,cluster_tool=None):
    with open(path_to_perfect_cluster,'rb') as f:
        pc=pickle.load(f)
    with open(path_to_candidate_cluster,'rb') as f:
        cc=pickle.load(f)
    if show_PC:
        if path_to_perfect_cluster:
            print(f'\nPerfect cluster: Total {len(pc)} results below\n')
            pprint.pprint(pc)
            print("\n----------\n")
        else:
            print("Please input the path to perfect_cluster.pkl!")
    if show_CC:
        if path_to_candidate_cluster:
            print(f'\nCandidate cluster: Total {len(cc)} results below\n')
            pprint.pprint(cc)
            print("\n----------\n")
        else:
            print("Please input the path to candidate_cluster.pkl!")
     
    def do_drep(dir_name,path_to_checkm_result):
        start=datetime.datetime.now()
        text=f'Starting dRep process at {dir_name}/ ...\n'
        for words in text:
            sys.stdout.write(words)
            sys.stdout.flush()
            time.sleep(0.05)
        
        if os.system('cd %s;\
                echo "genome,completeness,contamination,strain_heterogeneity" > sample_stat.csv;\
                sed "s/\t/,/g" %s > checkm;\
                for i in $(cat sample_hmq.tsv | awk -F "/" \'{print $NF}\' | sed "s/.fa//g");\
                do grep "$i," checkm \
                | awk -F "," \'{print $1".fa",$12,$13,$14}\' | sed "s/ /,/g";done >> sample_stat.csv;rm checkm' % (dir_name,path_to_checkm_result)):
            print("Generating sample_stat.csv failed! Checking script and try again!")
        else:
            with open(f'{dir_name}/drep.sh','w') as f:
                    f.write("#!/bin/bash\n\
                        dRep dereplicate \\\n\
                        ./ \\\n\
                        --processors 20 \\\n\
                        --length 50000 \\\n\
                        --completeness 75 \\\n\
                        --contamination 25 \\\n\
                        --genomes sample_hmq.tsv \\\n\
                        --genomeInfo sample_stat.csv \\\n\
                        --S_algorithm fastANI \\\n\
                        --P_ani 0.9 \\\n\
                        --S_ani 0.99 \\\n\
                        --cov_thresh 0.1 \\\n\
                        > drep.log 2>&1")
            if os.system(f'cd {dir_name};bash drep.sh'):
                print(f'Processing dRep failed! Checking {dir_name}/drep.log for more details!')
            else:
                print(f'Processing dRep finished! Checking {dir_name}/ for results!')
                end=datetime.datetime.now()
                running_time=end-start
                print("Running time: ",running_time)
                print("----------\n")
    
    if cluster_type and cluster_num:
        if path_to_perfect_cluster and cluster_type=='perfect':
            dir_name=f'{cluster_type}_cluster_{cluster_num}_ani95_drep'
            if os.system(f'ls {dir_name}'):
                os.system(f'mkdir {dir_name}')
                with open(f'{dir_name}/sample_hmq.tsv','w') as f:
                    for sample in pc[cluster_num]['members']:
                        f.writelines(f'{path_to_mag}/{sample}\n')
                do_drep(dir_name,path_to_checkm_result)
            else:
                print(f'Processing failed! The directory {dir_name} is already exist!')
        elif path_to_candidate_cluster and cluster_type=='candidate' and cluster_tool:
            if len(cc[cluster_num][cluster_tool])>1:
                dir_name=f'{cluster_type}_cluster_{cluster_num}_ani95_drep'
                if os.system(f'ls {dir_name}'):
                    os.system(f'mkdir {dir_name}')
                    with open(f'{dir_name}/sample_hmq.tsv','w') as f:
                        for sample in cc[cluster_num][cluster_tool]:
                            f.writelines(f'{path_to_mag}/{sample}\n')
                    do_drep(dir_name,path_to_checkm_result)
                else:
                    print(f'Processing failed! The directory {dir_name} is already exist!')
            else:
                print(f'Processing failed! The {cluster_tool} cluster of {cluster_type}\'s {cluster_num} is less than 2 members!')
        else:
            print("Please input correct path_to_perfect_cluster, cluster_type and cluster_tool(if choose candidate cluster) options!")
    
def main():
    chose_candidate_or_perfect_cluster(path_to_perfect_cluster='./perfect_cluster.pkl',path_to_candidate_cluster='./candidate_cluster.pkl',\
    path_to_checkm_result='/store/yhshao/META/metapi_result_cami/link_to_mag_vamb_no_graphbin/checkm_output/checkm.tsv',\
    path_to_mag='/store/yhshao/META/metapi_result_cami/link_to_mag_vamb_no_graphbin',\
    show_PC=True,show_CC=True,cluster_type='perfect',cluster_num='result9')#,cluster_tool='dRep')
    
if __name__=='__main__':
    main()

#Example output of this workflow:
# Perfect cluster: Total 9 results below

# {'result1': {'RG': 'CAMI_UT_2.metaspades.metabat2.bin.22.fa',
#              'members': ['CAMI_UT_2.metaspades.metabat2.bin.22.fa',
#                          'CAMI_UT_25.metaspades.concoct.bin.46_sub.fa',
#                          'CAMI_UT_3.metaspades.vamb.bin.13_sub.fa'],
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus '
#                          'helveticus']},
#  'result2': {'RG': 'CAMI_UT_24.metaspades.metabat2.bin.5.fa',
#              'members': ['CAMI_UT_22.metaspades.maxbin2.bin.2_sub.fa',
#                          'CAMI_UT_24.metaspades.metabat2.bin.5.fa',
#                          'CAMI_UT_25.metaspades.metabat2.bin.15.fa',
#                          'CAMI_UT_3.metaspades.metabat2.bin.6.fa'],
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus '
#                          'delbrueckii']},
#  'result3': {'RG': 'CAMI_UT_25.metaspades.metabat2.bin.7.fa',
#              'members': ['CAMI_UT_24.metaspades.metabat2.bin.15.fa',
#                          'CAMI_UT_25.metaspades.metabat2.bin.7.fa',
#                          'CAMI_UT_3.metaspades.maxbin2.bin.6_sub.fa',
#                          'CAMI_UT_5.metaspades.vamb.bin.8.fa'],
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Levilactobacillus;s__Levilactobacillus '
#                          'brevis']},
#  'result4': {'RG': 'CAMI_UT_3.metaspades.concoct.bin.25.fa',
#              'members': ['CAMI_UT_3.metaspades.concoct.bin.25.fa',
#                          'CAMI_UT_6.metaspades.concoct.bin.14.fa'],
#              'species': ['d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Bifidobacteriaceae;g__Bifidobacterium;s__Bifidobacterium '
#                          'vaginale_G']},
#  'result5': {'RG': 'CAMI_UT_24.metaspades.metabat2.bin.6_sub.fa',
#              'members': ['CAMI_UT_24.metaspades.metabat2.bin.6_sub.fa',
#                          'CAMI_UT_25.metaspades.maxbin2.bin.6.fa',
#                          'CAMI_UT_5.metaspades.metabat2.bin.6_sub.fa'],
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Latilactobacillus;s__Latilactobacillus '
#                          'curvatus']},
#  'result6': {'RG': 'CAMI_UT_2.metaspades.concoct.bin.21.fa',
#              'members': ['CAMI_UT_2.metaspades.concoct.bin.21.fa',
#                          'CAMI_UT_24.metaspades.metabat2.bin.10.fa',
#                          'CAMI_UT_25.metaspades.maxbin2.bin.1_sub.fa'],
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Limosilactobacillus;s__Limosilactobacillus '
#                          'fermentum']},
#  'result7': {'RG': 'CAMI_UT_2.metaspades.maxbin2.bin.1.fa',
#              'members': ['CAMI_UT_0.metaspades.concoct.bin.42.fa',
#                          'CAMI_UT_2.metaspades.maxbin2.bin.1.fa',
#                          'CAMI_UT_22.metaspades.maxbin2.bin.1.fa'],
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Ligilactobacillus;s__Ligilactobacillus '
#                          'salivarius']},
#  'result8': {'RG': 'CAMI_UT_5.metaspades.vamb.bin.14.fa',
#              'members': ['CAMI_UT_25.metaspades.vamb.bin.16.fa',
#                          'CAMI_UT_3.metaspades.vamb.bin.20.fa',
#                          'CAMI_UT_5.metaspades.vamb.bin.14.fa'],
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus '
#                          'amylovorus']},
#  'result9': {'RG': 'CAMI_UT_22.metaspades.concoct.bin.21.fa',
#              'members': ['CAMI_UT_2.metaspades.vamb.bin.2.fa',
#                          'CAMI_UT_22.metaspades.concoct.bin.21.fa'],
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lentilactobacillus;s__Lentilactobacillus '
#                          'buchneri']}}

# ----------

# Candidate cluster: Total 5 results below

# {'result1': {'RG': 'CAMI_UT_2.metaspades.maxbin2_graphbin_dastools.bin.18.fa',
#              'dRep': ['CAMI_UT_2.metaspades.maxbin2_graphbin_dastools.bin.18.fa',
#                       'CAMI_UT_22.metaspades.concoct_graphbin_dastools.bin.36_sub.fa'],
#              'galah': ['CAMI_UT_2.metaspades.maxbin2_graphbin_dastools.bin.18.fa'],
#              'same species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Secundilactobacillus;s__Secundilactobacillus '
#                               'paracollinoides']},
#  'result2': {'RG': 'CAMI_UT_21.metaspades.concoct_graphbin_dastools.bin.2.fa',
#              'dRep': ['CAMI_UT_0.metaspades.concoct_graphbin_dastools.bin.8_sub.fa',
#                       'CAMI_UT_2.metaspades.concoct_graphbin_dastools.bin.35_sub.fa',
#                       'CAMI_UT_21.metaspades.concoct_graphbin_dastools.bin.2.fa',
#                       'CAMI_UT_3.metaspades.metabat2_graphbin_dastools.bin.15_sub.fa'],
#              'galah': ['CAMI_UT_21.metaspades.concoct_graphbin_dastools.bin.2.fa'],
#              'same species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactiplantibacillus;s__Lactiplantibacillus '
#                               'plantarum']},
#  'result3': {'RG': 'CAMI_UT_22.metaspades.metabat2_graphbin_dastools.bin.10_sub.fa',
#              'dRep': ['CAMI_UT_22.metaspades.metabat2_graphbin_dastools.bin.10_sub.fa',
#                       'CAMI_UT_24.metaspades.maxbin2_graphbin_dastools.bin.8_sub.fa'],
#              'galah': ['CAMI_UT_22.metaspades.metabat2_graphbin_dastools.bin.10_sub.fa'],
#              'same species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lacticaseibacillus;s__Lacticaseibacillus '
#                               'paracasei']},
#  'result4': {'RG': 'CAMI_UT_3.metaspades.concoct_graphbin_dastools.bin.16.fa',
#              'dRep': ['CAMI_UT_25.metaspades.metabat2_graphbin_dastools.bin.18.fa',
#                       'CAMI_UT_3.metaspades.concoct_graphbin_dastools.bin.16.fa',
#                       'CAMI_UT_5.metaspades.maxbin2_graphbin_dastools.bin.11_sub.fa'],
#              'galah': ['CAMI_UT_3.metaspades.concoct_graphbin_dastools.bin.16.fa',
#                        'CAMI_UT_25.metaspades.metabat2_graphbin_dastools.bin.18.fa'],
#              'same species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Limosilactobacillus;s__Limosilactobacillus '
#                               'fermentum']},
#  'result5': {'RG': 'CAMI_UT_22.metaspades.maxbin2_graphbin_dastools.bin.2.fa',
#              'dRep': ['CAMI_UT_0.metaspades.concoct_graphbin_dastools.bin.26_sub.fa',
#                       'CAMI_UT_22.metaspades.maxbin2_graphbin_dastools.bin.2.fa'],
#              'galah': ['CAMI_UT_22.metaspades.maxbin2_graphbin_dastools.bin.2.fa'],
#              'same species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus '
#                               'delbrueckii']}}

# ----------

# Starting dRep process at perfect_cluster_result9_ani95_drep/ ...
# Processing dRep finished! Checking perfect_cluster_result9_ani95_drep/ for results!
# Running time:  0:00:09.826470
# ----------
