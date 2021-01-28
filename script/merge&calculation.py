import numpy as np
import pandas as pd
import glob

def cal_enrich(count_DF,column,input):
    count_DF[column+'_enrich'] = np.log10(((count_DF[column]+1)/(count_DF[column].sum(0)+143))/((count_DF[input]+1)/(count_DF[input].sum(0)+143)))
    return count_DF

def cal_mean(count_DF,mean,column1,column2):
    count_DF[mean] = (count_DF[column1]+count_DF[column2])/2
    return count_DF
def main():
#    nep_files = glob.glob('results/*/nep_mut.tsv')

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
    one_mut_sortdf = onemut_df.join(pos_df, lsuffix='_caller', rsuffix='_other').set_index(0).sort_index()

    muts_df = nep_df[nep_df["NEP_Mutation"].str.contains('-')]

    correlation = nep_df.corr(method='pearson')
    correlation.to_csv('results/nep_mut_correlation.csv')
    nep_df.to_csv('results/nep_mut.csv')
    one_mut_sortdf.to_csv('results/one_mut.csv')
    muts_df.to_csv('results/multi_mut.csv')




if __name__ == "__main__":
  main()