import os
import time
import datetime
import sys
import pprint

def choose_candidate_or_perfect_cluster(path_to_perfect_cluster=None,path_to_candidate_cluster=None,\
            show_PC=False,show_CC=False,cluster_type=None,cluster_num=None,cluster_tool=None):
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
     
    def do_drep(dir_name):
        start=datetime.datetime.now()
        text=f'Starting dRep process at {dir_name}/ ...\n'
        for words in text:
            sys.stdout.write(words)
            sys.stdout.flush()
            time.sleep(0.05)
        if os.system('cd %s;\
                echo "genome completeness contamination strain_heterogeneity N50 contig_number contig_length_sum contig_length_max" > sample_stat.csv;\
                for i in $(cat sample_hmq.tsv | cut -d "/" -f13);\
                do grep $i /hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/flye_test/bins_hmq.metaspades.dastools.gtdbtk.all.tsv \
                | awk \'{print $1,$(NF-11),$(NF-10),$(NF-9),$NF,$(NF-5),$(NF-4),$(NF-2)}\';done >> sample_stat.csv;\
                sed -i \'s/ /,/g\' sample_stat.csv' % dir_name):
            print("Generating sample_stat.csv failed! Checking script and try again!")
        else:
            with open(f'{dir_name}/drep.sh','w') as f:
                    f.write("#!/bin/bash\n\
                        dRep dereplicate \\\n\
                        ./ \\\n\
                        --processors 32 \\\n\
                        --length 50000 \\\n\
                        --completeness 75 \\\n\
                        --contamination 25 \\\n\
                        --genomes sample_hmq.tsv \\\n\
                        --genomeInfo sample_stat.csv \\\n\
                        --S_algorithm ANImf \\\n\
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
            os.system(f'mkdir {dir_name}')
            with open(f'{dir_name}/sample_hmq.tsv','w') as f:
                for sample in pc[cluster_num]['members']:
                    sample_id=sample.split('.')[0]
                    f.writelines(f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag/results/06.binning/bins/{sample_id}.metaspades.out/dastools/{sample}\n')
            do_drep(dir_name)
        elif path_to_candidate_cluster and cluster_type=='candidate' and cluster_tool:
            if len(cc[cluster_num][cluster_tool])>1:
                dir_name=f'{cluster_type}_cluster_{cluster_num}_ani95_drep'
                os.system(f'mkdir {dir_name}')
                with open(f'{dir_name}/sample_hmq.tsv','w') as f:
                    for sample in cc[cluster_num][cluster_tool]:
                        sample_id=sample.split('.')[0]
                        f.writelines(f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag/results/06.binning/bins/{sample_id}.metaspades.out/dastools/{sample}\n')
                do_drep(dir_name)
            else:
                print(f'Processing failed! The {cluster_tool} cluster of {cluster_type}\'s {cluster_num} is less than 2 members!')
        else:
            print("Please input correct path_to_perfect_cluster, cluster_type and cluster_tool(if choose candidate cluster) options!")
    
def main():
    choose_candidate_or_perfect_cluster(path_to_perfect_cluster='./perfect_cluster.pkl',path_to_candidate_cluster='./candidate_cluster.pkl',\
            show_PC=False,show_CC=True,cluster_type='candidate',cluster_num='result16',cluster_tool='galah')
    
if __name__=='__main__':
    main()
