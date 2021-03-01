import pandas as pd 
import os, sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pprint

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

def venn_plot_of_RG(drep_RG,galah_RG):
    venn2([set(drep_RG),set(galah_RG)],set_labels=('drep_RG','galah_RG'))

def get_drep_or_galah_only_cluster_dict(drep_cluster_table,drep_RG,galah_cluster_table,galah_RG,show_drep=False,show_galah=False):
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
                col1=table.columns[0]
                col2=table.columns[1]
            else:
                col1=table.columns[1]
                col2=table.columns[0]
            for i in diff:
                start_index=list(table[col1]).index(i)
                cluster_members_index=[j for j in range(len(table)) if table[col2][start_index]==table[col2][j]]
                cluster[i]=list(table[col1][cluster_members_index[0]:(cluster_members_index[-1]+1)])
            return cluster

    if show_drep:
        print('dRep:')
        pprint.pprint(mkdict(drep_RG,galah_RG,drep_cluster_table))
    if show_galah:
        print('galah:')
        pprint.pprint(mkdict(galah_RG,drep_RG,galah_cluster_table))

def get_intersection_cluster_dict_of_drep_or_galah(drep_cluster_table,drep_RG,galah_cluster_table,galah_RG,show_inter_cluster_drep=False,show_inter_cluster_galah=False,show_diff_of_inter=False):
    def mkdict_drep(RG1,RG2,table):
        diff=[i for i in RG1 if i not in RG2]
        col1=table.columns[0]
        col2=table.columns[1]
        inter_cluster_drep={}.fromkeys(diff)
        for i in RG1:
            if i not in diff:
                start_index=list(table[col1]).index(i)
                cluster_members_index=[j for j in range(len(table)) if table[col2][start_index]==table[col2][j]]
                inter_cluster_drep[i]=list(table[col1][cluster_members_index[0]:(cluster_members_index[-1]+1)])
        return inter_cluster_drep
    
    def mkdict_galah(RG1,RG2,table):
        diff=[i for i in RG1 if i not in RG2]
        col1=table.columns[0]
        col2=table.columns[1]
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
    
    if show_inter_cluster_drep:
        print('dRep:')
        pprint.pprint(mkdict_drep(drep_RG,galah_RG,drep_cluster_table))
    if show_inter_cluster_galah:
        print('galah:')
        pprint.pprint(mkdict_galah(galah_RG,drep_RG,galah_cluster_table))
    if show_diff_of_inter:
        inter=[i for i in drep_RG if i in galah_RG]
        inter_cluster_drep=mkdict_drep(drep_RG,galah_RG,drep_cluster_table)
        inter_cluster_galah=mkdict_galah(galah_RG,drep_RG,galah_cluster_table)
        pprint.pprint([('*',k,'<dRep>:',inter_cluster_drep[k],'<galah>:',inter_cluster_galah[k]) for k in inter if sorted(inter_cluster_drep[k]) != sorted(inter_cluster_galah[k])])
        print('**************')

    
def main():
    assembler_name=['megahit','metaspades']
    binner_name=['dastools','metabat2']
    plt.ion()
    for i in assembler_name:
        for j in binner_name:
            path_to_drep_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/metapiin_temp/results/10.dereplicate/mags/hmq.bins.{i}.{j}.drep.out/data_tables/Cdb.csv'
            path_to_galah_cluster_table=f'/hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/metapiin_temp/results/10.dereplicate_galah/{i}_{j}/output-cluster-definition.tsv'
            path_to_drep_RG_dir=f'/hwfssz1/ST_META/P18Z10200N0127_MA/shaoyanhan/metapitest1/metapiin_temp/results/10.dereplicate/mags/hmq.bins.{i}.{j}.drep.out/dereplicated_genomes'
            a1,a2=read_in_drep_cluster_table(path_to_drep_cluster_table,path_to_drep_RG_dir)
            b1,b2=read_in_galah_cluster_table(path_to_galah_cluster_table)
#             print(a1,a2,b1,b2)
            print(f'{i}_{j}')
            venn_plot_of_RG(a2,b2)
            plt.pause(.2)
            print('**************')
            get_drep_or_galah_only_cluster_dict(a1,a2,b1,b2,show_drep=True,show_galah=True)
            print('**************')
            get_intersection_cluster_dict_of_drep_or_galah(a1,a2,b1,b2,show_inter_cluster_drep=False,show_inter_cluster_galah=False,show_diff_of_inter=True)
    plt.ioff()
    plt.show()
if __name__=='__main__':
    main()
