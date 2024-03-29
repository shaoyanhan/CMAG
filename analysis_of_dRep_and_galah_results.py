import pandas as pd 
import os, sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pprint
import pickle

def read_in_drep_cluster_table(path_to_drep_cluster_table,path_to_drep_RG_dir):
    drep_RG=[]
    drep_RG=os.listdir(path_to_drep_RG_dir)
    if '.ipynb_checkpoints' in drep_RG:
        drep_RG.remove('.ipynb_checkpoints')
    return pd.read_csv(path_to_drep_cluster_table,sep=','),drep_RG

def read_in_galah_cluster_table(path_to_galah_cluster_table):
    galah_cluster_table=[]
    galah_cluster_table=pd.read_csv(path_to_galah_cluster_table,sep='\t',header=None,names=['RG','member'])
    for i in range(len(galah_cluster_table)):
        galah_cluster_table['RG'][i]=galah_cluster_table['RG'][i].split('/')[-1]
        galah_cluster_table['member'][i]=galah_cluster_table['member'][i].split('/')[-1]
    galah_RG=[]
    [galah_RG.append(i) for i in galah_cluster_table['RG'] if i not in galah_RG]
    return galah_cluster_table,galah_RG

def read_in_taxon_table(path_to_taxon_table):
    return pd.read_csv(path_to_taxon_table,sep='\t')

def venn_plot_of_RG(drep_RG,galah_RG):
    venn2([set(drep_RG),set(galah_RG)],set_labels=('drep_RG','galah_RG'))

def get_drep_or_galah_only_cluster_dict(drep_cluster_table,drep_RG,galah_cluster_table,galah_RG,\
            show_drep=False,show_galah=False,return_drep=False,return_galah=False,return_same_mamber_cluster=False):
    def mkdict(RG1,RG2,table):
        table_type=len(table.columns)
        diff=[i for i in RG1 if i not in RG2]    
        if len(diff)==0:
            if table_type>2:
                print('There are no drep only clusters, because drep_RG is a subset of galah_RG!')
            else:
                print('There are no galah only clusters, because galah_RG is a subset of drep_RG!')
        else:   
            cluster={}.fromkeys(diff)
            if table_type>2:
                col1='genome' #table.columns[0]
                col2='secondary_cluster' #table.columns[1]
            else:
                col1='member' #table.columns[1]
                col2='RG' #table.columns[0]
            for i in diff:
                start_index=list(table[col1]).index(i)
                cluster_members_index=[j for j in range(len(table)) if table[col2][start_index]==table[col2][j]]
                cluster[i]=list(table[col1][cluster_members_index[0]:(cluster_members_index[-1]+1)])
            return cluster
    
    drep_result=mkdict(drep_RG,galah_RG,drep_cluster_table)
    galah_result=mkdict(galah_RG,drep_RG,galah_cluster_table)
    same_member_cluster={}
    for i in drep_result:
        for j in galah_result:
            if sorted(drep_result[i])==sorted(galah_result[j]):
                same_member_cluster[i]=drep_result[i]
                break
    if same_member_cluster:
        print('\n**************')
        print(f'There are {len(same_member_cluster)} clusters between {len(drep_result)} drep_only and {len(galah_result)} galah_only clusters with different RG but same members!\n')
        pprint.pprint(same_member_cluster)
        print('**************\n')
    else:
        print('\n**************')
        print(f'There are no clusters with different RG but same members!')
        print('**************\n')
    
    if return_same_mamber_cluster:
        return same_member_cluster
    
    returns=[]
    if show_drep:
        print('\n**************')
        print(f'dRep_only: total= {len(drep_result)} results')
        pprint.pprint(drep_result)
        print('**************\n')
    if show_galah:
        print('\n**************')
        print(f'galah_only: total= {len(galah_result)} results')
        pprint.pprint(galah_result)
        print('**************\n')
    if return_drep:
        returns.append(drep_result)
    if return_galah:
        returns.append(galah_result)
    if returns:
        return returns

def get_intersection_cluster_dict_of_drep_or_galah(drep_cluster_table,drep_RG,galah_cluster_table,galah_RG,\
    show_inter_cluster_drep=False,show_inter_cluster_galah=False,show_diff_of_inter=False,return_drep=False,return_galah=False,return_diff_or_same=None):
    def mkdict_drep(RG1,RG2,table):
        inter=[i for i in RG1 if i in RG2]
        col1='genome' #table.columns[0]
        col2='secondary_cluster' #table.columns[1]
        inter_cluster_drep={}.fromkeys(inter)
        for i in inter:
            start_index=list(table[col1]).index(i)
            cluster_members_index=[j for j in range(len(table)) if table[col2][start_index]==table[col2][j]]
            inter_cluster_drep[i]=list(table[col1][cluster_members_index[0]:(cluster_members_index[-1]+1)])
        return inter_cluster_drep
    
    def mkdict_galah(RG1,RG2,table):
        diff=[i for i in RG1 if i not in RG2]
        col1='RG'#table.columns[0]
        col2='member'#table.columns[1]
        i=0
        j=0
        depth_of_cluster=[[]]
        while i<(len(table)-1):
            if table[col1][i]!=table[col1][i+1]:
                depth_of_cluster.append([i-j,j+1])
                i=i+1
                j=0
            else:
                i=i+1
                j=j+1
        if j==0:
            depth_of_cluster.append([i,j+1])
        else:
            depth_of_cluster.append([i-j,j+1])
        depth_of_cluster.pop(0)
        index_to_delet=[]
        for i in diff:
            index_to_delet.append(RG1.index(i))
        depth_of_cluster=[depth_of_cluster[i] for i in range(len(depth_of_cluster)) if i not in index_to_delet]
        inter_cluster_galah={}
        for i in depth_of_cluster:
            inter_cluster_galah[table[col1][i[0]]]=list(table[col2][i[0]:(i[0]+i[1])])
        return inter_cluster_galah
    
    returns=[]
    if show_inter_cluster_drep:
        result=mkdict_drep(drep_RG,galah_RG,drep_cluster_table)
        print('\n**************')
        print(f'dRep_inter: total= {len(result)} results')
        pprint.pprint(result)
        print('**************\n')
    if show_inter_cluster_galah:
        result=mkdict_galah(galah_RG,drep_RG,galah_cluster_table)
        print('\n**************')
        print(f'galah_inter: total= {len(result)} results')
        pprint.pprint(result)
        print('**************\n')
    if show_diff_of_inter:
        inter=[i for i in drep_RG if i in galah_RG]
        inter_cluster_drep=mkdict_drep(drep_RG,galah_RG,drep_cluster_table)
        inter_cluster_galah=mkdict_galah(galah_RG,drep_RG,galah_cluster_table)
        result=[('*',k,'<dRep>:',inter_cluster_drep[k],'<galah>:',inter_cluster_galah[k]) for k in inter if sorted(inter_cluster_drep[k]) != sorted(inter_cluster_galah[k])]
        print('\n**************')
        print(f'difference between clusters of intersections: total= {len(result)} results')
        pprint.pprint(result)
        print('**************\n')     
    if return_drep:
        returns.append(mkdict_drep(drep_RG,galah_RG,drep_cluster_table))
    if return_galah:
        returns.append(mkdict_galah(galah_RG,drep_RG,galah_cluster_table))
    if return_diff_or_same:
        inter=[i for i in drep_RG if i in galah_RG]
        inter_cluster_drep=mkdict_drep(drep_RG,galah_RG,drep_cluster_table)
        inter_cluster_galah=mkdict_galah(galah_RG,drep_RG,galah_cluster_table)
        if return_diff_or_same=='diff':
            return [(k,inter_cluster_drep[k],inter_cluster_galah[k]) for k in inter if sorted(inter_cluster_drep[k]) != sorted(inter_cluster_galah[k])]
        elif return_diff_or_same=='same':
            dict_list=[{k:inter_cluster_drep[k]} for k in inter if sorted(inter_cluster_drep[k]) == sorted(inter_cluster_galah[k])]
            result=dict_list.pop()
            for i in range(len(dict_list)):
                result.update(dict_list.pop())
            same_member_cluster=get_drep_or_galah_only_cluster_dict(drep_cluster_table,drep_RG,galah_cluster_table,galah_RG,return_same_mamber_cluster=True)
            if same_member_cluster:
                result.update(same_member_cluster)
                print('\n**************')
                print(f'\nAdding {len(same_member_cluster)} clusters with different RG but same members to the "same" result set!\n')
                print('**************\n')
            return result
        else:
            print('Please input the selection (\'diff\'/\'same\') you want to return!')
    if returns:
        return returns

    
def taxonomy_verification(taxon_table,diff=None,same=None,show_diff=None,show_same=None,check_cluster_replication=False,write_to_file=None):
    def taxon_diff(taxon_table,diff):
        same_taxon={}
        diff_taxon={}
        for i in range(len(diff)):
            members=diff[i][1] if len(diff[i][1])>len(diff[i][2]) else diff[i][2]
            taxon=[]
            for j in members:
                taxon.append(taxon_table['GTDB classification'][list(taxon_table['user_genome']).index(j)])
            if len(list(set(taxon)))==1:
                result={}
                result['RG']=diff[i][0]
                result['same species']=list(set(taxon))
                result['dRep']=diff[i][1]
                result['galah']=diff[i][2]
                same_taxon[f'result{len(same_taxon)+1}']=result
            else:
                result={}
                result['RG']=diff[i][0]
                result['each\'s species']=taxon
                result['dRep']=diff[i][1]
                result['galah']=diff[i][2]
                diff_taxon[f'result{len(diff_taxon)+1}']=result
        return same_taxon,diff_taxon

    def taxon_same(taxon_table,same):
        same_taxon_single={}
        same_taxon_multi={}
        diff_taxon={}
        single_RG_count=0
        for i in same.keys():
            if len(same[i])==1:
                taxon=[]
                taxon.append(taxon_table['GTDB classification'][list(taxon_table['user_genome']).index(i)])
                result={}
                result['RG']=i
                result['species']=taxon
                result['member']=i
                same_taxon_single[f'result{len(same_taxon_single)+1}']=result
                single_RG_count=single_RG_count+1
            else:
                taxon=[]
                for j in same[i]:
                    taxon.append(taxon_table['GTDB classification'][list(taxon_table['user_genome']).index(j)])
                if len(list(set(taxon)))==1:
                    result={}
                    result['RG']=i
                    result['species']=list(set(taxon))
                    result['members']=same[i]
                    same_taxon_multi[f'result{len(same_taxon_multi)+1}']=result
                else:
                    result={}
                    result['RG']=i
                    result['each\'s species']=taxon
                    result['dRep']=same[i]
                    diff_taxon[f'result{len(diff_taxon)+1}']=result
        return same_taxon_single,same_taxon_multi,single_RG_count,diff_taxon
    
    if show_diff:
        same_taxon,diff_taxon=taxon_diff(taxon_table,diff)
        print('\n**************')
        print('Show: RG same but members not same results')
        print('-------------')
        print(f'| S |  There are total {len(same_taxon)+len(diff_taxon)} results below:')
        print(f'| T |  Cluster\'s members from same species: {len(same_taxon)}; ')
        print(f'| A |  Cluster\'s members from different species: {len(diff_taxon)}; ')
        print('| T |')
        if show_diff=='all':
            print('-------------')
            print(f'<same species>: total= {len(same_taxon)} results')
            pprint.pprint(same_taxon)
            print('-------------')
            print(f'<different species>: total= {len(diff_taxon)} results')
            pprint.pprint(diff_taxon)
        print('**************\n')
        
    if show_same:
        same_taxon_single,same_taxon_multi,single_RG_count,diff_taxon=taxon_same(taxon_table,same)
        print('\n**************')
        print('Show: RG and members same results')
        print('-------------')
        print(f'| S |  There are total {single_RG_count+len(same_taxon_multi)+len(diff_taxon)} results below:')
        print(f'| T |  Cluster with only one member: {single_RG_count}; ')
        print(f'| A |  Cluster with at least two members and are from same species: {len(same_taxon_multi)}; ')
        print(f'| T |  Cluster\'s members from different species: {len(diff_taxon)}; ')
        print('-------------')
        if show_same=='simple':
            print(f'Omitting {single_RG_count} results of one member clusters')
            print('-------------')
            print(f'<same species>: total= {len(same_taxon_multi)} results')
            pprint.pprint(same_taxon_multi)
            print('-------------')
            print(f'<different species>: total= {len(diff_taxon)} results')
            pprint.pprint(diff_taxon)
        elif show_same=='all':
            print(f'<one member>: total= {single_RG_count} results')
            pprint.pprint(same_taxon_single)
            print('-------------')
            print(f'<same species>: total= {len(same_taxon_multi)} results')
            pprint.pprint(same_taxon_multi)
            print('-------------')
            print(f'<different species>: total= {len(diff_taxon)} results')
            pprint.pprint(diff_taxon)         
        else:
            print('Please select \'simple\' or \'all\' mode to print results')
        print('**************\n')
            
    def species_point_to_result_dict(same_taxon):
        species_point_to_result={}
        replication_count=0
        for i in same_taxon:
            species_name=same_taxon[i]['species'][0]
            if species_name in species_point_to_result:
                species_point_to_result[species_name].append(i)
                replication_count=replication_count+1
            else:
                species_point_to_result[species_name]=[i]
        return species_point_to_result,replication_count
            
    def replication_result_dict(same_taxon,species_point_to_result):
        replication_result={}
        for i in species_point_to_result:
            per_species_result=species_point_to_result[i]
            if len(per_species_result)>1:
                result={}
                result['result count']=len(per_species_result)
                result['result number of <one member> or <same species>']=per_species_result
                result['RG']=[]
                for j in per_species_result:
                    result['RG'].append(same_taxon[j]['RG'])
                result['species']=i
                replication_result[f'result{len(replication_result)+1}']=result
        return replication_result
    
    if check_cluster_replication:
        same_taxon_single,same_taxon_multi=taxon_same(taxon_table,same)[0:2]
        species_point_to_result_single,replication_count_single=species_point_to_result_dict(same_taxon_single)
        species_point_to_result_multi,replication_count_multi=species_point_to_result_dict(same_taxon_multi)
#         print(species_point_to_result_single,replication_count_single)
#         print('-------------')
#         print(species_point_to_result_multi,replication_count_multi)
        if replication_count_single>0:
            replication_result_single=replication_result_dict(same_taxon_single,species_point_to_result_single)
        if replication_count_multi>0:
            replication_result_multi=replication_result_dict(same_taxon_multi,species_point_to_result_multi)
        print('-------------')
        print('<check results>')
        if replication_count_single>0 and replication_count_multi==0:
            print(f'There are total {len(replication_result_single)} replicate species between the set of <one member> clusters!')
            print(f'There are no replicate species between the set of <same species> clusters!')
            pprint.pprint(replication_result_single)
        elif replication_count_single==0 and replication_count_multi>0:
            print(f'There are no replicate species between the set of <one member> clusters!')
            print(f'There are total {len(replication_result_multi)} replicate species between the set of <same species> clusters!')
            pprint.pprint(replication_result_multi)
        elif replication_count_single>0 and replication_count_multi>0:
            print(f'There are total {len(replication_result_single)} and {len(replication_result_multi)} replicate species between the set of <one member> and <same species> clusters!')
            print(f'<one member>: {len(replication_result_single)} results')
            pprint.pprint(replication_result_single)
            print('-------------')
            print(f'<same species>: {len(replication_result_multi)} results')
            pprint.pprint(replication_result_multi)
        else:
            print('There are no replicate species between all same species clusters.')
        
    if write_to_file:
        if write_to_file=='CC_and_PC':
            same_taxon=taxon_diff(taxon_table,diff)[0]
            same_taxon_multi=taxon_same(taxon_table,same)[1]
            with open('candidate_cluster.pkl','wb') as f:
                pickle.dump(same_taxon,f)
            with open('perfect_cluster.pkl','wb') as f:
                pickle.dump(same_taxon_multi,f)
            print('Candidate and perfect cluster were already written to files: candidate_cluster.pkl perfect_cluster.pkl!')
        else:
            print("Please input option to write files!")
        
def main():
#     assembler_name=['megahit','metaspades']
#     binner_name=['dastools','metabat2']
#     plt.ion()
#     for i in assembler_name:
#         for j in binner_name:
#             path_to_drep_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/metapiin_temp/results/10.dereplicate/mags/hmq.bins.{i}.{j}.drep.out/data_tables/Cdb.csv'
#             path_to_galah_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/metapiin_temp/results/10.dereplicate_galah/{i}_{j}/output-cluster-definition.tsv'
#             path_to_drep_RG_dir=f'/hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/metapiin_temp/results/10.dereplicate/mags/hmq.bins.{i}.{j}.drep.out/dereplicated_genomes'
#             path_to_taxon_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/metapiin_temp/results/09.classify/report/bins_hmq.{i}.{j}.gtdbtk.all.tsv'
#             a1,a2=read_in_drep_cluster_table(path_to_drep_cluster_table,path_to_drep_RG_dir)
#             b1,b2=read_in_galah_cluster_table(path_to_galah_cluster_table)
#             c=read_in_taxon_table(path_to_taxon_table)
# #             print(a1,a2,b1,b2)
#             print(f'{i}_{j}')
#             venn_plot_of_RG(a2,b2)
#             plt.pause(.2)
# #             print('**************')
# #             get_drep_or_galah_only_cluster_dict(a1,a2,b1,b2,show_drep=True,show_galah=True)
#             print('**************')
#             get_intersection_cluster_dict_of_drep_or_galah(a1,a2,b1,b2,show_inter_cluster_drep=False,show_inter_cluster_galah=False,show_diff_of_inter=True)
# #             print('**************')
# #             same=get_intersection_cluster_dict_of_drep_or_galah(a1,a2,b1,b2,return_diff_or_same='same')
# #             diff=get_intersection_cluster_dict_of_drep_or_galah(a1,a2,b1,b2,return_diff_or_same='diff')
# #             taxonomy_verification(c,diff,same,check_cluster_replication=True)
#             print('**************\n')
#     plt.ioff()
#     plt.show()

#CAMI
    path_to_drep_cluster_table=f'/store/yhshao/META/metapi_result_cami/drep_vamb_no_graphbin_beyond_50_completeness_no_contamination/data_tables/Cdb.csv'
    path_to_galah_cluster_table=f'/store/yhshao/META/metapi_result_cami/galah_vamb_no_graphbin_beyond_50_completeness_no_contamination/cluster.tsv'
    path_to_drep_RG_dir=f'/store/yhshao/META/metapi_result_cami/drep_vamb_no_graphbin_beyond_50_completeness_no_contamination/dereplicated_genomes/'
    path_to_taxon_table=f'/store/yhshao/META/metapi_result_cami/gtdbtk.bac120.summary.tsv'

#ANI95
#     path_to_drep_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools_drep/ANI95/data_tables/Cdb.csv'
#     path_to_galah_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools.tsv'
#     path_to_drep_RG_dir=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools_drep/ANI95/dereplicated_genomes'
#     path_to_taxon_table=f'/zfssz3/ST_META/ST_META_CD/PROJECT/P18Z10200N0127_ZJ/vagina/vagina_mag/results/09.classify/report/bins_hmq.metaspades.dastools.gtdbtk.all.tsv'

#as_mt_qh
#     path_to_drep_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools_as_mt_qh_drep_galah/drep/data_tables/Cdb.csv'
#     path_to_galah_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools_as_mt_qh_drep_galah/galah/cluster-definition.tsv'
#     path_to_drep_RG_dir=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools_as_mt_qh_drep_galah/drep/dereplicated_genomes'
#     path_to_taxon_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools_as_mt_qh_drep_galah/cluster_mag_dastools_as_mt_qh.taxonmy.tsv'

# ANI99
#     path_to_drep_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools_drep/ANI99/data_tables/Cdb.csv'
#     path_to_galah_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools_drep/galah_ani99/cluster-definition.tsv'
#     path_to_drep_RG_dir=f'/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/vagina/vagina_mag_analysis/assay/cluster/cluster_mag_dastools_drep/ANI99/dereplicated_genomes'
#     path_to_taxon_table=f'/zfssz3/ST_META/ST_META_CD/PROJECT/P18Z10200N0127_ZJ/vagina/vagina_mag/results/09.classify/report/bins_hmq.metaspades.dastools.gtdbtk.all.tsv'


    a1,a2=read_in_drep_cluster_table(path_to_drep_cluster_table,path_to_drep_RG_dir)
    b1,b2=read_in_galah_cluster_table(path_to_galah_cluster_table)
    c=read_in_taxon_table(path_to_taxon_table)
    plt.ion()
    venn_plot_of_RG(a2,b2)
    plt.ioff()
    plt.show()
    
    get_drep_or_galah_only_cluster_dict(a1,a2,b1,b2,show_drep=True,show_galah=True)
    get_intersection_cluster_dict_of_drep_or_galah(a1,a2,b1,b2,show_inter_cluster_drep=True,show_inter_cluster_galah=True,show_diff_of_inter=True)
    same=get_intersection_cluster_dict_of_drep_or_galah(a1,a2,b1,b2,return_diff_or_same='same')
    diff=get_intersection_cluster_dict_of_drep_or_galah(a1,a2,b1,b2,return_diff_or_same='diff')
    taxonomy_verification(c,diff,same,show_diff='all',show_same='all',check_cluster_replication=True,write_to_file='CC_and_PC')
    
if __name__=='__main__':
    main()
    
    
#Example output of this workflow:
#Draw the venn picture of dRep and galah RG 

# **************
# There are 3 clusters between 3 drep_only and 3 galah_only clusters with different RG but same members!

# {'CAMI_UT_2.metaspades.maxbin2.bin.1.fa': ['CAMI_UT_0.metaspades.concoct.bin.42.fa',
#                                            'CAMI_UT_2.metaspades.maxbin2.bin.1.fa',
#                                            'CAMI_UT_22.metaspades.maxbin2.bin.1.fa'],
#  'CAMI_UT_22.metaspades.concoct.bin.21.fa': ['CAMI_UT_2.metaspades.vamb.bin.2.fa',
#                                              'CAMI_UT_22.metaspades.concoct.bin.21.fa'],
#  'CAMI_UT_5.metaspades.vamb.bin.14.fa': ['CAMI_UT_25.metaspades.vamb.bin.16.fa',
#                                          'CAMI_UT_3.metaspades.vamb.bin.20.fa',
#                                          'CAMI_UT_5.metaspades.vamb.bin.14.fa']}
# **************


# **************
# dRep_only: total= 3 results
# {'CAMI_UT_2.metaspades.maxbin2.bin.1.fa': ['CAMI_UT_0.metaspades.concoct.bin.42.fa',
#                                            'CAMI_UT_2.metaspades.maxbin2.bin.1.fa',
#                                            'CAMI_UT_22.metaspades.maxbin2.bin.1.fa'],
#  'CAMI_UT_22.metaspades.concoct.bin.21.fa': ['CAMI_UT_2.metaspades.vamb.bin.2.fa',
#                                              'CAMI_UT_22.metaspades.concoct.bin.21.fa'],
#  'CAMI_UT_5.metaspades.vamb.bin.14.fa': ['CAMI_UT_25.metaspades.vamb.bin.16.fa',
#                                          'CAMI_UT_3.metaspades.vamb.bin.20.fa',
#                                          'CAMI_UT_5.metaspades.vamb.bin.14.fa']}
# **************


# **************
# galah_only: total= 3 results
# {'CAMI_UT_2.metaspades.vamb.bin.2.fa': ['CAMI_UT_2.metaspades.vamb.bin.2.fa',
#                                         'CAMI_UT_22.metaspades.concoct.bin.21.fa'],
#  'CAMI_UT_22.metaspades.maxbin2.bin.1.fa': ['CAMI_UT_22.metaspades.maxbin2.bin.1.fa',
#                                             'CAMI_UT_0.metaspades.concoct.bin.42.fa',
#                                             'CAMI_UT_2.metaspades.maxbin2.bin.1.fa'],
#  'CAMI_UT_3.metaspades.vamb.bin.20.fa': ['CAMI_UT_3.metaspades.vamb.bin.20.fa',
#                                          'CAMI_UT_5.metaspades.vamb.bin.14.fa',
#                                          'CAMI_UT_25.metaspades.vamb.bin.16.fa']}
# **************


# **************
# dRep_inter: total= 14 results
# {'CAMI_UT_2.metaspades.concoct.bin.21.fa': ['CAMI_UT_2.metaspades.concoct.bin.21.fa',
#                                             'CAMI_UT_24.metaspades.metabat2.bin.10.fa',
#                                             'CAMI_UT_25.metaspades.maxbin2.bin.1_sub.fa'],
#  'CAMI_UT_2.metaspades.maxbin2.bin.14.fa': ['CAMI_UT_2.metaspades.maxbin2.bin.14.fa'],
#  'CAMI_UT_2.metaspades.metabat2.bin.10_sub.fa': ['CAMI_UT_2.metaspades.metabat2.bin.10_sub.fa'],
#  'CAMI_UT_2.metaspades.metabat2.bin.22.fa': ['CAMI_UT_2.metaspades.metabat2.bin.22.fa',
#                                              'CAMI_UT_25.metaspades.concoct.bin.46_sub.fa',
#                                              'CAMI_UT_3.metaspades.vamb.bin.13_sub.fa'],
#  'CAMI_UT_2.metaspades.vamb.bin.17.fa': ['CAMI_UT_2.metaspades.vamb.bin.17.fa'],
#  'CAMI_UT_24.metaspades.metabat2.bin.5.fa': ['CAMI_UT_22.metaspades.maxbin2.bin.2_sub.fa',
#                                              'CAMI_UT_24.metaspades.metabat2.bin.5.fa',
#                                              'CAMI_UT_25.metaspades.metabat2.bin.15.fa',
#                                              'CAMI_UT_3.metaspades.metabat2.bin.6.fa'],
#  'CAMI_UT_24.metaspades.metabat2.bin.6_sub.fa': ['CAMI_UT_24.metaspades.metabat2.bin.6_sub.fa',
#                                                  'CAMI_UT_25.metaspades.maxbin2.bin.6.fa',
#                                                  'CAMI_UT_5.metaspades.metabat2.bin.6_sub.fa'],
#  'CAMI_UT_25.metaspades.metabat2.bin.7.fa': ['CAMI_UT_24.metaspades.metabat2.bin.15.fa',
#                                              'CAMI_UT_25.metaspades.metabat2.bin.7.fa',
#                                              'CAMI_UT_3.metaspades.maxbin2.bin.6_sub.fa',
#                                              'CAMI_UT_5.metaspades.vamb.bin.8.fa'],
#  'CAMI_UT_3.metaspades.concoct.bin.25.fa': ['CAMI_UT_3.metaspades.concoct.bin.25.fa',
#                                             'CAMI_UT_6.metaspades.concoct.bin.14.fa'],
#  'CAMI_UT_3.metaspades.concoct.bin.3.fa': ['CAMI_UT_3.metaspades.concoct.bin.3.fa'],
#  'CAMI_UT_6.metaspades.concoct.bin.12.fa': ['CAMI_UT_6.metaspades.concoct.bin.12.fa'],
#  'CAMI_UT_6.metaspades.concoct.bin.2.fa': ['CAMI_UT_6.metaspades.concoct.bin.2.fa'],
#  'CAMI_UT_6.metaspades.maxbin2.bin.1.fa': ['CAMI_UT_6.metaspades.maxbin2.bin.1.fa'],
#  'CAMI_UT_6.metaspades.maxbin2.bin.3.fa': ['CAMI_UT_6.metaspades.maxbin2.bin.3.fa']}
# **************


# **************
# galah_inter: total= 14 results
# {'CAMI_UT_2.metaspades.concoct.bin.21.fa': ['CAMI_UT_2.metaspades.concoct.bin.21.fa',
#                                             'CAMI_UT_24.metaspades.metabat2.bin.10.fa',
#                                             'CAMI_UT_25.metaspades.maxbin2.bin.1_sub.fa'],
#  'CAMI_UT_2.metaspades.maxbin2.bin.14.fa': ['CAMI_UT_2.metaspades.maxbin2.bin.14.fa'],
#  'CAMI_UT_2.metaspades.metabat2.bin.10_sub.fa': ['CAMI_UT_2.metaspades.metabat2.bin.10_sub.fa'],
#  'CAMI_UT_2.metaspades.metabat2.bin.22.fa': ['CAMI_UT_2.metaspades.metabat2.bin.22.fa',
#                                              'CAMI_UT_3.metaspades.vamb.bin.13_sub.fa',
#                                              'CAMI_UT_25.metaspades.concoct.bin.46_sub.fa'],
#  'CAMI_UT_2.metaspades.vamb.bin.17.fa': ['CAMI_UT_2.metaspades.vamb.bin.17.fa'],
#  'CAMI_UT_24.metaspades.metabat2.bin.5.fa': ['CAMI_UT_24.metaspades.metabat2.bin.5.fa',
#                                              'CAMI_UT_3.metaspades.metabat2.bin.6.fa',
#                                              'CAMI_UT_25.metaspades.metabat2.bin.15.fa',
#                                              'CAMI_UT_22.metaspades.maxbin2.bin.2_sub.fa'],
#  'CAMI_UT_24.metaspades.metabat2.bin.6_sub.fa': ['CAMI_UT_24.metaspades.metabat2.bin.6_sub.fa',
#                                                  'CAMI_UT_5.metaspades.metabat2.bin.6_sub.fa',
#                                                  'CAMI_UT_25.metaspades.maxbin2.bin.6.fa'],
#  'CAMI_UT_25.metaspades.metabat2.bin.7.fa': ['CAMI_UT_25.metaspades.metabat2.bin.7.fa',
#                                              'CAMI_UT_24.metaspades.metabat2.bin.15.fa',
#                                              'CAMI_UT_5.metaspades.vamb.bin.8.fa',
#                                              'CAMI_UT_3.metaspades.maxbin2.bin.6_sub.fa'],
#  'CAMI_UT_3.metaspades.concoct.bin.25.fa': ['CAMI_UT_3.metaspades.concoct.bin.25.fa',
#                                             'CAMI_UT_6.metaspades.concoct.bin.14.fa'],
#  'CAMI_UT_3.metaspades.concoct.bin.3.fa': ['CAMI_UT_3.metaspades.concoct.bin.3.fa'],
#  'CAMI_UT_6.metaspades.concoct.bin.12.fa': ['CAMI_UT_6.metaspades.concoct.bin.12.fa'],
#  'CAMI_UT_6.metaspades.concoct.bin.2.fa': ['CAMI_UT_6.metaspades.concoct.bin.2.fa'],
#  'CAMI_UT_6.metaspades.maxbin2.bin.1.fa': ['CAMI_UT_6.metaspades.maxbin2.bin.1.fa'],
#  'CAMI_UT_6.metaspades.maxbin2.bin.3.fa': ['CAMI_UT_6.metaspades.maxbin2.bin.3.fa']}
# **************


# **************
# difference between clusters of intersections: total= 0 results
# []
# **************


# **************
# There are 3 clusters between 3 drep_only and 3 galah_only clusters with different RG but same members!

# {'CAMI_UT_2.metaspades.maxbin2.bin.1.fa': ['CAMI_UT_0.metaspades.concoct.bin.42.fa',
#                                            'CAMI_UT_2.metaspades.maxbin2.bin.1.fa',
#                                            'CAMI_UT_22.metaspades.maxbin2.bin.1.fa'],
#  'CAMI_UT_22.metaspades.concoct.bin.21.fa': ['CAMI_UT_2.metaspades.vamb.bin.2.fa',
#                                              'CAMI_UT_22.metaspades.concoct.bin.21.fa'],
#  'CAMI_UT_5.metaspades.vamb.bin.14.fa': ['CAMI_UT_25.metaspades.vamb.bin.16.fa',
#                                          'CAMI_UT_3.metaspades.vamb.bin.20.fa',
#                                          'CAMI_UT_5.metaspades.vamb.bin.14.fa']}
# **************


# **************

# Adding 3 clusters with different RG but same members to the "same" result set!

# **************


# **************
# Show: RG same but members not same results
# -------------
# | S |  There are total 0 results below:
# | T |  Cluster's members from same species: 0; 
# | A |  Cluster's members from different species: 0; 
# | T |
# -------------
# <same species>: total= 0 results
# {}
# -------------
# <different species>: total= 0 results
# {}
# **************


# **************
# Show: RG and members same results
# -------------
# | S |  There are total 17 results below:
# | T |  Cluster with only one member: 8; 
# | A |  Cluster with at least two members and are from same species: 9; 
# | T |  Cluster's members from different species: 0; 
# -------------
# <one member>: total= 8 results
# {'result1': {'RG': 'CAMI_UT_6.metaspades.maxbin2.bin.3.fa',
#              'member': 'CAMI_UT_6.metaspades.maxbin2.bin.3.fa',
#              'species': ['d__Bacteria;p__Firmicutes_C;c__Negativicutes;o__Veillonellales;f__Dialisteraceae;g__Dialister_A;s__Dialister_A '
#                          'pneumosintes']},
#  'result2': {'RG': 'CAMI_UT_6.metaspades.concoct.bin.2.fa',
#              'member': 'CAMI_UT_6.metaspades.concoct.bin.2.fa',
#              'species': ['d__Bacteria;p__Actinobacteriota;c__Coriobacteriia;o__Coriobacteriales;f__Eggerthellaceae;g__Denitrobacterium;s__Denitrobacterium '
#                          'detoxificans']},
#  'result3': {'RG': 'CAMI_UT_6.metaspades.concoct.bin.12.fa',
#              'member': 'CAMI_UT_6.metaspades.concoct.bin.12.fa',
#              'species': ['d__Bacteria;p__Actinobacteriota;c__Coriobacteriia;o__Coriobacteriales;f__Eggerthellaceae;g__Slackia;s__Slackia '
#                          'heliotrinireducens']},
#  'result4': {'RG': 'CAMI_UT_3.metaspades.concoct.bin.3.fa',
#              'member': 'CAMI_UT_3.metaspades.concoct.bin.3.fa',
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Mycoplasmatales;f__Mycoplasmoidaceae;g__Ureaplasma;s__Ureaplasma '
#                          'parvum']},
#  'result5': {'RG': 'CAMI_UT_6.metaspades.maxbin2.bin.1.fa',
#              'member': 'CAMI_UT_6.metaspades.maxbin2.bin.1.fa',
#              'species': ['d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella '
#                          'intermedia']},
#  'result6': {'RG': 'CAMI_UT_2.metaspades.vamb.bin.17.fa',
#              'member': 'CAMI_UT_2.metaspades.vamb.bin.17.fa',
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Limosilactobacillus;s__Limosilactobacillus '
#                          'mucosae']},
#  'result7': {'RG': 'CAMI_UT_2.metaspades.metabat2.bin.10_sub.fa',
#              'member': 'CAMI_UT_2.metaspades.metabat2.bin.10_sub.fa',
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lacticaseibacillus;s__Lacticaseibacillus '
#                          'paracasei']},
#  'result8': {'RG': 'CAMI_UT_2.metaspades.maxbin2.bin.14.fa',
#              'member': 'CAMI_UT_2.metaspades.maxbin2.bin.14.fa',
#              'species': ['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lacticaseibacillus;s__Lacticaseibacillus '
#                          'rhamnosus']}}
# -------------
# <same species>: total= 9 results
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
# -------------
# <different species>: total= 0 results
# {}
# **************

# -------------
# <check results>
# There are no replicate species between all same species clusters.
# Candidate and perfect cluster were already written to files: candidate_cluster.pkl perfect_cluster.pkl!
