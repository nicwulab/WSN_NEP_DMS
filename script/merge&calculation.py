import numpy as np
import pandas as pd
import glob

def cal_RF(count_DF,column_name,input):
    count_DF[column_name+'_RF'] = np.log10(((count_DF[column_name]+1)/(count_DF[column_name].sum(0)+len(count_DF[column_name])))/((count_DF[input]+1)/(count_DF[input].sum(0)+len(count_DF[column_name]))))
    return count_DF

def cal_sum(count_DF,sum_name,column1,column2):
    count_DF[sum_name] = (count_DF[column1]+count_DF[column2])
    return count_DF

def cal_mean(count_DF,mean,column1,column2):
    count_DF[mean] = (count_DF[column1]+count_DF[column2])/2
    return count_DF
def main():
#merge the sample nep_mutation and output the single mutation file(one_mut.csv) and multi-mutation file(multi_mut.csv)

    nep1_df = pd.read_csv('results/Sample4_TGACCAAT/nep_mut.tsv', header=0, sep= '\t')
    nep2_df = pd.read_csv('results/Sample5_ACAGTGAT/nep_mut.tsv', header=0, sep= '\t')
    nep3_df = pd.read_csv('results/Sample6_GCCAATAT/nep_mut.tsv', header=0, sep= '\t')
    nep4_df = pd.read_csv('results/Sample7_CAGATCAT/nep_mut.tsv', header=0, sep= '\t')
    nep_df1 = pd.merge(nep1_df,nep2_df,how='outer',on='NEP_Mutation')
    nep_df1.columns = ['NEP_Mutation', 'input1', 'input2']
    nep_df2 = pd.merge(nep3_df,nep4_df,how='outer',on='NEP_Mutation')
    nep_df2.columns = ['NEP_Mutation', 'output1', 'output2']
    nep_df = pd.merge(nep_df1,nep_df2,how='outer',on='NEP_Mutation')
    onemut_df = nep_df[~nep_df["NEP_Mutation"].str.contains('-')]
    pos_df = onemut_df.NEP_Mutation.str.extract('(\d+)')
    one_mut_sortdf = onemut_df.join(pos_df, lsuffix='_caller', rsuffix='_other')\
        .set_index(0)\
        .fillna(0)\
        .sort_index()
    muts_df = nep_df[nep_df["NEP_Mutation"].str.contains('-')] #filter by string contain '-'
    correlation = one_mut_sortdf.corr(method='pearson')
    correlation.to_csv('results/nep_mut_correlation.csv')

    Analyzed_onemut_df = cal_sum(one_mut_sortdf,'input','input1','input2')
    cal_RF(Analyzed_onemut_df,'output1','input')
    cal_RF(Analyzed_onemut_df, 'output2', 'input')
    cal_mean(Analyzed_onemut_df,'average_RF','output1_RF','output2_RF')
    Analyzed_onemut_df.to_csv('results/Analyzed_onemut.csv')
#    nep_df.to_csv('results/nep_mut.csv')
#    one_mut_sortdf.to_csv('results/one_mut.csv')
#    muts_df.to_csv('results/multi_mut.csv')




if __name__ == "__main__":
  main()