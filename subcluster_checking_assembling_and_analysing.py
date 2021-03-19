import os
import time
import datetime
import sys
import pprint

def get_cluster_dict(path_to_drep_cluster_table):
    with open(path_to_drep_cluster_table,'r') as f:
        cluster_table=f.readlines()
    for i in range(len(cluster_table)):
        cluster_table[i]=cluster_table[i].rstrip('\n')
    cluster_dict={}
    for i in cluster_table[1:]:
        j=i.split(',')
        if j[1] not in cluster_dict:
            cluster_dict[j[1]]=[j[0]]
        else:
            cluster_dict[j[1]].append(j[0])
    print(f'Removing {len([j for j in cluster_dict if len(cluster_dict[j])<2])} <one member> clusters')
    [cluster_dict.pop(i) for i in [j for j in cluster_dict if len(cluster_dict[j])<2]]
    print(f'Total left {len(cluster_dict)} clusters')
    return cluster_dict

def check_subcluster_animf99(cluster_dict):
    failed_count=0
    for cluster in cluster_dict:
        os.system("mkdir drep")

        with open('./drep/sample_hmq.tsv','w') as f:
            for sample in cluster_dict[cluster]:
                sample_id=sample.split('.')[0]
                f.writelines(f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag/results/06.binning/bins/{sample_id}.metaspades.out/dastools/{sample}\n')

        os.system("echo \"genome completeness contamination strain_heterogeneity N50 contig_number contig_length_sum contig_length_max\" > ./drep/sample_stat.csv;\
            for i in $(cat ./drep/sample_hmq.tsv | cut -d \"/\" -f13);\
            do grep $i /hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/flye_test/bins_hmq.metaspades.dastools.gtdbtk.all.tsv \
            | awk '{print $1,$6,$7,$8,$17,$12,$13,$15}';done >> ./drep/sample_stat.csv;\
            sed -i 's/ /,/g' ./drep/sample_stat.csv")

        with open('./drep/drep.sh','w') as f:
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

        os.system("chmod +x ./drep/drep.sh")

        text=f'Processing dRep of cluster {cluster} ...\n'
        for words in text:
            sys.stdout.write(words)
            sys.stdout.flush()
            time.sleep(0.3)

        start=datetime.datetime.now()
        while os.system("cd drep/;bash drep.sh")==0:
            print(f'Cluster {cluster} finished dRep process!')
            break

        f=os.popen("ls ./drep/dereplicated_genomes/ | wc -l")
        if f.read()[0]=='1':
            print(f'Cluster {cluster} passed examination!')
            failed=False
        else:
            print(f'There are multiple subsets in cluster {cluster}')
            failed=True

        if failed:
            with open('./failed_cluster.log','a') as f:
                f.write(f'\n<<<<<<<<<<Cluster {cluster}>>>>>>>>>>\n')
            os.system("cat ./drep/data_tables/Cdb.csv >> ./failed_cluster.log")
            failed_count=failed_count+1

        with open('./all_drep.log','a') as f:
            f.write(f'<<<<<<<<<<Cluster {cluster}>>>>>>>>>>\n')
        os.system("cat ./drep/sample_stat.csv ./drep/drep.log >> ./all_drep.log;\
                    ls ./drep/dereplicated_genomes/ >>./all_drep.log;\
                    rm -rf ./drep/")

        end=datetime.datetime.now()
        running_time=end-start
        with open('./all_drep.log','a') as f:
            f.write(f'Running time: {running_time}\n\n')
        print("Running time: ",running_time)
        print("----------\n")

    if failed_count:
        print(f'All checking process finished, {failed_count} clusters failed. Checking drep.log and failed_cluster.log for more details.')
    else:
        print("All checking process finished, all clusters passed. Checking drep.log for more details.")

def make_RG_statistical_file(path_to_drep_RG_dir,cluster_dict):
    drep_RG=[]
    drep_RG=os.listdir(path_to_drep_RG_dir)
    if '.ipynb_checkpoints' in drep_RG:
        drep_RG.remove('.ipynb_checkpoints')
    # print(len(drep_RG))
    RG_to_cluster=[]
    for cluster in cluster_dict:
        for sample in cluster_dict[cluster]:
            if sample in drep_RG:
                RG_to_cluster.append([sample,cluster])
                drep_RG.remove(sample)
                break
    with open('RG_to_cluster.tsv','w') as f:
        f.write('RG Cluster \n')
        for i in RG_to_cluster:
            f.write(' '.join(i)+' \n')
            
    text=f'Writting RG\'s statistical file to RG_stat.csv ...\n'
    for words in text:
        sys.stdout.write(words)
        sys.stdout.flush()
        time.sleep(0.05)
    if os.system("echo \"completeness contamination strain_heterogeneity N50 contig_number contig_length_sum contig_length_max\" > RG_stat.tsv;\
                for i in $(cat RG_to_cluster.tsv | cut -d \" \" -f1 | sed -n '2,$p');\
                do grep $i /hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/flye_test/bins_hmq.metaspades.dastools.gtdbtk.all.tsv \
                | awk '{print $6,$7,$8,$17,$12,$13,$15}';done >> RG_stat.tsv;\
                paste RG_to_cluster.tsv RG_stat.tsv > RG_stat.csv;\
                sed -i 's/ /,/g' RG_stat.csv;\
                rm RG_stat.tsv RG_to_cluster.tsv"):
        print("Writting failed! Please checking your script and running again!")
    else:
        print("Writting finished! Checking ./RG_stat.csv for more details!")

def subcluster_flye_assembling_and_checkm_checking(cluster_dict):
    os.system("mkdir -p ./flye/assemblies/")
    os.system("mkdir -p ./flye/logs/")
    with open('./flye/assembly_stat.csv','a') as f:
        f.write("Genome,Total_length,Fragments,N50,Largest_frg,Scaffolds,Mean_coverage\n")
    succeed=True
    for cluster in cluster_dict:
        os.system("mkdir ./flye_temp/")
        with open('./flye_temp/flye.sh','a') as f:
            f.write("#!/bin/bash\nflye --subassemblies \\\n")
            for sample in cluster_dict[cluster]:
                sample_id=sample.split('.')[0]
                f.writelines(f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag/results/06.binning/bins/{sample_id}.metaspades.out/dastools/{sample} \\\n')
            f.write("--threads 32 \\\n-o ./")
            
        start=datetime.datetime.now()
        text=f'Assembling cluster {cluster}\'s MAGs using flye ...\n'
        for words in text:
            sys.stdout.write(words)
            sys.stdout.flush()
            time.sleep(0.05)
        if os.system("cd ./flye_temp/;chmod +x flye.sh;bash flye.sh>out.log 2>&1"):
            print(f'Cluster {cluster}\'s flye subassemblies processing failed! Checking ./flye_temp/out.log and ./flye_temp/flye.log for more details.')
            succeed=False
            break
        else:
            print(f'Cluster {cluster}\'s flye subassemblies processing successed!')
            os.system(f'cd ./flye_temp/;\
                        tail -8 flye.log | head -6 | cut -d ":" -f2 | xargs echo \"{cluster}.fasta\" >>../flye/assembly_stat.csv;\
                        mv assembly.fasta ../flye/assemblies/{cluster}.fasta;\
                        cat flye.log >../flye/logs/{cluster}.log')
        os.system("rm -rf ./flye_temp/")
        
        end=datetime.datetime.now()
        running_time=end-start
        print("Running time: ",running_time)
        print("----------\n")
    os.system("sed -i 's/ /,/g' ./flye/assembly_stat.csv")
    print('All flye processes finished! Checking ./flye/ or ./flye_temp/(if failed) for assemblies, statistical table and logs.')
    
    if succeed:
        start=datetime.datetime.now()
        text=f'Starting checkm process ...\n'
        for words in text:
            sys.stdout.write(words)
            sys.stdout.flush()
            time.sleep(0.2)
        if os.system("cd ./flye/assemblies/;checkm lineage_wf --tab_table --file ../checkm_assemblies.tsv -x fasta -t 32 ./ ../checkm_assemblies_output"):
            print('CheckM processing failed! Checking ./flye/checkm_assemblies_output/checkm.log for more details')
        else:
            print('CheckM processing finished! Checking ./flye/checkm_assemblies.tsv for results.')
            end=datetime.datetime.now()
            running_time=end-start
            print("Running time: ",running_time)
            print("----------\n")
    else:
        print('Subassemblies processing failed, automatically shutting down the checkM process!')
    print('All flye and checkm processes finished! Checking ./flye for more details.')
        
def main():
    path_to_drep_cluster_table='/hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/candidate_cluster_result6_ani95_drep/data_tables/Cdb.csv'
    path_to_drep_RG_dir='/hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/candidate_cluster_result6_ani95_drep/dereplicated_genomes/'
    cluster_dict=get_cluster_dict(path_to_drep_cluster_table)
    pprint.pprint(cluster_dict)
    check_subcluster_animf99(cluster_dict)
    subcluster_flye_assembling_and_checkm_checking(cluster_dict)
    make_RG_statistical_file(path_to_drep_RG_dir,cluster_dict)
    os.system("mv RG_stat.csv ./flye/")
if __name__=='__main__':
    main()
    
        
