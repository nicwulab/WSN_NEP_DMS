import numpy as np
import pandas as pd
import glob

def cal_enri1(count_DF,column_name,input):
    count_DF[column_name+'_enri'] = ((count_DF[column_name]+1)/(99653+2072))/((count_DF[input]+1)/(80095+2072))
    count_DF[column_name + '_RF'] = np.log10(count_DF[column_name+'_enri'])
    wt_DF = count_DF[count_DF.NEP_Mutation=='WT'][column_name + '_RF'].iloc[0]
    count_DF[column_name + '_RF'] = (count_DF[column_name + '_RF'] - wt_DF)
    return count_DF
def cal_enri2(count_DF,column_name,input):
    count_DF[column_name+'_enri'] = ((count_DF[column_name]+1)/(83632+2072))/((count_DF[input]+1)/(91389+2072))
    count_DF[column_name + '_RF'] = np.log10(count_DF[column_name+'_enri'])
    wt_DF = count_DF[count_DF.NEP_Mutation=='WT'][column_name + '_RF'].iloc[0]
    count_DF[column_name + '_RF'] = (count_DF[column_name + '_RF'] - wt_DF)
    return count_DF
def cal_min(count_DF,min,column1,column2):
    count_DF[min] = count_DF[[column1, column2]].min(axis=1)
    return count_DF
def main():
#merge the sample nep_mutation and output the single mutation file(one_mut.csv) and multi-mutation file(multi_mut.csv)

    nep1_df = pd.read_csv('results/Sample4_TGACCAAT/nep_mut_ns1wt.tsv', header=0, sep= '\t')
    nep2_df = pd.read_csv('results/Sample5_ACAGTGAT/nep_mut_ns1wt.tsv', header=0, sep= '\t')
    nep3_df = pd.read_csv('results/Sample6_GCCAATAT/nep_mut_ns1wt.tsv', header=0, sep= '\t')
    nep4_df = pd.read_csv('results/Sample7_CAGATCAT/nep_mut_ns1wt.tsv', header=0, sep= '\t')
    nep_df1 = pd.merge(nep1_df,nep2_df,how='outer',on='NEP_Mutation')
    nep_df1.columns = ['NEP_Mutation', 'input1', 'input2']
    nep_df2 = pd.merge(nep3_df,nep4_df,how='outer',on='NEP_Mutation')
    nep_df2.columns = ['NEP_Mutation', 'output1', 'output2']
    nep_df = pd.merge(nep_df1,nep_df2,how='outer',on='NEP_Mutation')
    count_df = nep_df.loc[(nep_df['input1']>10) & (nep_df['input2']>10)]

    pos_df = count_df.NEP_Mutation.str.extract('(\d+)')
    one_mut_sortdf = count_df.join(pos_df, lsuffix='_caller', rsuffix='_other')\
        .set_index(0)\
        .fillna(0)\
        .sort_index()
    cal_enri1(one_mut_sortdf,'output1','input1')
    cal_enri2(one_mut_sortdf,'output2','input2')
    cal_min(one_mut_sortdf,'average_RF','output1_RF','output2_RF')
    one_mut_sortdf.to_csv('results/Analyzed_onemut_ns1wt.csv')
if __name__ == "__main__":
  main()


