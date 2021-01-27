import numpy as np
import pandas as pd

def cal_enrich(count_DF,column,input):
    count_DF[column+'_enrich'] = np.log10(((count_DF[column]+1)/(count_DF[column].sum(0)+143))/((count_DF[input]+1)/(count_DF[input].sum(0)+143)))
    return count_DF

def cal_mean(count_DF,mean,column1,column2):
    count_DF[mean] = (count_DF[column1]+count_DF[column2])/2
    return count_DF
def main():
    nep1_df = pd.read_csv('script/nep1_mut.tsv', header=0, sep= '\t')
    nep2_df = pd.read_csv('script/nep2_mut.tsv', header=0, sep= '\t')
    nep3_df = pd.read_csv('script/nep3_mut.tsv', header=0, sep= '\t')
    nep4_df = pd.read_csv('script/nep4_mut.tsv', header=0, sep= '\t')
    nep_df1 = pd.merge(nep1_df,nep2_df,how='inner',on='NEP_Mutation')
    nep_df1.columns = ['NEP_Mutation', 'input1', 'input2']
    nep_df2 = pd.merge(nep3_df,nep4_df,how='inner',on='NEP_Mutation')
    nep_df2.columns = ['NEP_Mutation', 'output1', 'output2']
    nep_df = pd.merge(nep_df1,nep_df2,how='inner',on='NEP_Mutation')
    cal_mean(nep_df,'average_input','input1','input2')
    correlation = nep_df.corr(method='pearson')
    correlation.to_csv('result/nep_mut_correlation.csv')
    nep_df.to_csv('result/nep_mut.csv')




if __name__ == "__main__":
  main()