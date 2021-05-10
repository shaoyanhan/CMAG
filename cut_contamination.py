import os
import time
import datetime
import sys
import pprint
import pickle

def get_mag_multicopy_info_dict(mag_list,gene_copy_type,path_to_checkm_ext_file,show=False):
    print(f'\n***Writing mag\'s multicopy gene accession code to mag_multicopy_info_dictionary, mag list: {mag_list} ***\n')
    mag_multicopy_info={}
    for i in range(len(mag_list)):
        mag_multicopy_info[mag_list[i]]={}
        for j in range(len(gene_copy_type)):
            f=os.popen("grep -w \"%s\" %s | awk -F \"'%s': \" '{print $2}' | cut -d '[' -f2 | cut -d ']' -f1 | cut -d ',' -f1-" % (mag_list[i],path_to_checkm_ext_file,gene_copy_type[j]))#"'PF02686', 'PF01807'\n"
            accessions=f.read()
            mag_multicopy_info[mag_list[i]][gene_copy_type[j]]=accessions.lstrip('\"\'').rstrip('\'\n').split('\', \'')#{'CAMI_UT_6.metaspades.maxbin2.bin.2':{'GCN2': ['PF02686', 'PF01807'],'GCN3': [''],'GCN4': [''],'GCN5+': ['']}}
    if show:
        print("mag_multicopy_info:\n")
        pprint.pprint(mag_multicopy_info)
    return mag_multicopy_info



def read_in_mag_file(mag_name,mag_suffix,path_to_mag_folder):
    print('\n***Reading in mag file ...***\n')
    with open(f'{path_to_mag_folder}/{mag_name}.{mag_suffix}','r') as f:
        mag_file=f.readlines()#['>NODE_27_length_288659_cov_40.979365\n','CGGGCTACGCATGTGTCCGTTTCCTTGGCCATGGAGAGTGCGCCATTAAA\n',...]
    contigs={}#{'>NODE_27_length_288659_cov_40.979365':'CGGGCTACGCATGTGTCCGTTTCCTTGGCCATGGAGAGT...',...}
    i=0
    while i<len(mag_file):
        if mag_file[i][0]=='>':
            contig_name=mag_file[i].rstrip()[1:]
            contigs[contig_name]=mag_file[i+1].rstrip()
            i=i+2
        else:
            contigs[contig_name]=contigs[contig_name]+mag_file[i].rstrip()
            i=i+1
    return contigs

def append_contig_to_seq_interval_dict(cut_seq_info,hmm_result):#{'NODE_697_length_11598_cov_4.457252': [[4881, 6023], [1962, 3236], [1962, 3236]],...}
    hmm_contig_name=hmm_result.pop(0)#'NODE_697_length_11598_cov_4.457252_7'
    strain_type=hmm_result.pop()#'-1'
    contig_length=int(hmm_contig_name.split('_')[3])#11598
    contig_name='_'.join(hmm_contig_name.split('_')[0:6])#'NODE_697_length_11598_cov_4.457252'
    seq_interval=[]
#     if strain_type=='-1':#不管是正负链，checkm的结果都是显示正链的区间
#         seq_interval.append([contig_length-int(hmm_result[1])+1,contig_length-int(hmm_result[0])+1])#5576,6718 => [11598-6718+1,11598-5576+1]
#     else:
    seq_interval.append([int(hmm_result[0]),int(hmm_result[1])])
    if contig_name in cut_seq_info:
        cut_seq_info[contig_name].extend(seq_interval)
    else:
        cut_seq_info[contig_name]=seq_interval
    return cut_seq_info

def find_min_max(internal1,internal2):#[4,6],[3,5] => [3,6]
    internal1.extend(internal2)
    internal1.sort()
    return [internal1[0],internal1[3]]

def merge_overlap(internal_array):#[[4,7],[1,2],[3,5],[1,3],[8,10]] => [[1,7],[8,10]]
    merged_interval=[]
    while(internal_array):
        internal1=internal_array.pop(0)
        while True:
            i=0
            change=0
            while(i<len(internal_array)):
                if internal1[1]<internal_array[i][0]:
                    i=i+1
                else:
                    if internal1[0]>internal_array[i][1]:
                        i=i+1
                    else:
                        internal1=find_min_max(internal1,internal_array[i])
                        change=change+1
                        internal_array.pop(i)
            if change==0:
                break
        merged_interval.append(internal1)
    for i in range(len(merged_interval)-1):#[[8,10],[4,5],[1,2]] => [[1,2],[4,5],[8,10]]
        for j in range(len(merged_interval)-i-1):
            if merged_interval[j][0]>merged_interval[j+1][0]:
                temp=merged_interval[j]
                merged_interval[j]=merged_interval[j+1]
                merged_interval[j+1]=temp
    return merged_interval
        
def do_cutting(raw_contigs,cut_seq_info):# raw_contigs: {'NODE_1_length_78_cov_1': 'TTAACAGCGCG...',...}, cut_seq_info: {'NODE_1_length_78_cov_1': [[1,2],[3,4]],...}
    def length_check(raw_contigs,contig_name,raw_contigs_seq,head,end,split_count):
        if end-head+1>=12:# smallest gene length: start codon, stop codon and two amino acid
            split_count=split_count+1
            raw_contigs[f'{contig_name}_cutted_{split_count}']=raw_contigs_seq[(head-1):end]
        return raw_contigs,split_count
    
    print("\n***Starting cutting process ...***\n")
    for contig_name in cut_seq_info:# 'NODE_1_length_78_cov_1'
        raw_contigs_seq=raw_contigs.pop(contig_name)# 'TTAA...'
        contig_length=int(contig_name.split('_')[3])# 78
        interval_list=cut_seq_info[contig_name]# [[1,2],[3,4]]
        split_count=0
        
        if len(interval_list)==1:
            if interval_list[0][0]>1:
                raw_contigs,split_count=length_check(raw_contigs,contig_name,raw_contigs_seq,1,interval_list[0][0]-1,split_count)
            if interval_list[0][1]<contig_length:
                raw_contigs,split_count=length_check(raw_contigs,contig_name,raw_contigs_seq,interval_list[0][1]+1,contig_length,split_count)
                
        else:
            for i in range(len(interval_list)):
                if i==0:# the first interval [1,2]
                    if interval_list[i][0]==1:# [1,2]
                        continue
                    else:# the first interval is like [2,3]
                        raw_contigs,split_count=length_check(raw_contigs,contig_name,raw_contigs_seq,1,interval_list[i][0]-1,split_count)
                        continue
                elif i==len(interval_list)-1 and interval_list[i][1]!=contig_length:# the last interval [3,4]
                        raw_contigs,split_count=length_check(raw_contigs,contig_name,raw_contigs_seq,interval_list[i][1]+1,contig_length,split_count)
                raw_contigs,split_count=length_check(raw_contigs,contig_name,raw_contigs_seq,interval_list[i-1][1]+1,interval_list[i][0]-1,split_count)
                    
        if split_count:
            print(f'The contig {contig_name} was split to {split_count} segments. Interval list of this contig\'s contaminations: {interval_list}\n')
        else:
            print('\n----------\n')
            print(f'The contig {contig_name} was highly contaminated and thrown away! Interval list of this contig\'s contaminations: {interval_list}\n')
            print('\n----------\n')
    print("\n***Cutting process finished!***\n")
    return raw_contigs
    
def write_in_cutted_mag(cutted_contigs,mag_name,path_to_save_the_cutted_mag,mag_suffix):
    if os.system(f'ls {path_to_save_the_cutted_mag}'):
        os.system(f'mkdir -p {path_to_save_the_cutted_mag}')
    path_to_cutted_file=f'{path_to_save_the_cutted_mag}/{mag_name}_cutted.{mag_suffix}'
    with open(path_to_cutted_file,'w') as f:
        for contig_name in cutted_contigs:
            f.write('>'+contig_name+'\n')
            base_count=0
            for bases in cutted_contigs[contig_name]:
                if base_count==60:
                    f.write('\n')
                    base_count=0
                f.write(bases)
                base_count=base_count+1
            f.write('\n')
    print(f'Cutted contigs were written to file: {path_to_cutted_file} !')
    
def main():
    mag_suffix='fa'
    mag_list=['CAMI_UT_0.metaspades.metabat2.graphbin2.bin.18']#CAMI_UT_6.metaspades.maxbin2.bin.2']
    path_to_mag_folder='/store/yhshao/META/metapi_result_cami/link_to_mag_vamb_graphbin2'
    path_to_checkm_bins_folder='/store/yhshao/META/metapi_result_cami/checkm_output/bins'
    path_to_checkm_storage_folder='/store/yhshao/META/metapi_result_cami/checkm_output/storage'
    path_to_checkm_ext_file='/store/yhshao/META/metapi_result_cami/checkm_output/storage/bin_stats_ext.tsv'
    path_to_save_the_cutted_mag='/store/yhshao/META/cutted_mag'
    gene_copy_type=['GCN2','GCN3','GCN4','GCN5+']
    all_cut_mode=False
    all_cut_threshold='1e-10'
   
    mag_multicopy_info=get_mag_multicopy_info_dict(mag_list,gene_copy_type,path_to_checkm_ext_file,show=True)
    for mag_name in mag_multicopy_info:#'CAMI_UT_6.metaspades.maxbin2.bin.2',...

        start=datetime.datetime.now()
        text=f'\n<<< Starting contamination cutting process with {mag_name} ... >>>\n'
        for words in text:
            sys.stdout.write(words)
            sys.stdout.flush()
            time.sleep(0.01)

        cut_seq_info={}#{'NODE_697_length_11598_cov_4.457252': [4881, 6023, 1962, 3236, 1962, 3236],...}
        print("\n***Filling contigs and contaminations seq information into the cut_seq_info dictionary ...***\n")
        for multicopy_type in mag_multicopy_info[mag_name]:#'GCN2','GCN3','GCN4','GCN5+'
            for domain_accession in mag_multicopy_info[mag_name][multicopy_type]:#'PF02686', 'PF01807',...
                if domain_accession:
                    if all_cut_mode:
                        if all_cut_threshold:
                            # 'NODE_697_length_11598_cov_4.457252_7 5576 6718 -1' E-value: $7=>"1.5e-113" , cut all of the genes that E-value below the all_cut_threshold
                            f=os.popen("grep -w \"%s\" %s/%s/hmmer.analyze.txt | awk '{if($7<%s)print $1,$24,$26,$28}' | sed -n '2,$p'" % (domain_accession,path_to_checkm_bins_folder,mag_name,all_cut_threshold))
                        else:
                            print("Please input the all_cut_threshold!")
                    else:
                        if multicopy_type=='GCN5+':
                            print(f'GCN5+ exists in {mag_name}, just choosing TOP 2~5 contamination to cut!')
                        # 'NODE_697_length_11598_cov_4.457252_7 5576 6718 -1', Don't care about the E-value, only cut the called multicopy, e.g. GCN2 type only cut the second multicopy.
                        f=os.popen("grep -w \"%s\" %s/%s/hmmer.analyze.txt | awk '{print $1,$24,$26,$28}' | sed -n '2,%sp'" % (domain_accession,path_to_checkm_bins_folder,mag_name,multicopy_type[3]))

                    hmm_result_list=f.read().split('\n')
                    #'NODE_697_length_11598_cov_4.457252_7 5576 6718 -1\n' => ['NODE_697_length_11598_cov_4.457252_7 5576 6718 -1','']
                    #'NODE_11329_length_918_cov_1.783591_3 844 918 -1\nNODE_436_length_18290_cov_4.347993_1 1 75 1\n' => ['NODE_11329_length_918_cov_1.783591_3 844 918 -1','NODE_436_length_18290_cov_4.347993_1 1 75 1','']
                    hmm_result_list.pop()
                    for result in hmm_result_list:
                        hmm_result=result.split(' ')#'NODE_697_length_11598_cov_4.457252_7 5576 6718 -1' => ['NODE_697_length_11598_cov_4.457252_7','5576','6718','-1']
                        cut_seq_info=append_contig_to_seq_interval_dict(cut_seq_info,hmm_result)

        print("\n***Filling process finished!***\n")

        print("\n***Starting overlap merging process ...***\n")
        for contig_interval_list in cut_seq_info:#{'NODE_1177_length_6645_cov_37.049635': [[4656, 5495]],...}
            if len(cut_seq_info[contig_interval_list])>1:
                cut_seq_info[contig_interval_list]=merge_overlap(cut_seq_info[contig_interval_list])
        print("\n***Merging process finished!***\n")

        raw_contigs=read_in_mag_file(mag_name,mag_suffix,path_to_mag_folder)
        cutted_contigs=do_cutting(raw_contigs,cut_seq_info)
        write_in_cutted_mag(cutted_contigs,mag_name,path_to_save_the_cutted_mag,mag_suffix)

        print(f'Contamination cutting with {mag_name} finished!\n')
        end=datetime.datetime.now()
        running_time=end-start
        print("Running time: ",running_time)
        print("----------\n")

    print("All processes finished!")
    
if __name__=='__main__':
    main()
