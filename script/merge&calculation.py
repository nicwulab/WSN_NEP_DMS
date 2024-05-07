import numpy as np
import pandas as pd
import glob

def cal_enri(count_DF,column_name,input):
    count_DF[column_name+'_enri'] = ((count_DF[column_name]+1)/(count_DF[column_name].sum(0)+len(count_DF[column_name])))/((count_DF[input]+1)/(count_DF[input].sum(0)+len(count_DF[column_name])))
    count_DF[column_name + '_RF'] = np.log10(count_DF[column_name+'_enri'])
    wt_DF = count_DF[count_DF.NEP_Mutation=='WT'][column_name + '_RF'].iloc[0]
    count_DF[column_name + '_RF'] = (count_DF[column_name + '_RF'] - wt_DF)
    return count_DF

def cal_sum(count_DF,sum_name,column1,column2):
    count_DF[sum_name] = (count_DF[column1]+count_DF[column2])
    return count_DF

def cal_mean(count_DF,mean,column1,column2):
    count_DF[mean] = (count_DF[column1]+count_DF[column2])/2
    return count_DF

def cal_min(count_DF,min,column1,column2):
    count_DF[min] = count_DF[[column1, column2]].min(axis=1)
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

    count_df = nep_df.loc[(nep_df['input1']>10) & (nep_df['input2']>10)]

    count_df['Mutation_num'] = count_df['NEP_Mutation'].apply(lambda x: x.count('-')+1)
    count_df.to_csv('results/NEP_count_table.csv')

    onemut_df = nep_df[~nep_df["NEP_Mutation"].str.contains('-')]
    pos_df = onemut_df.NEP_Mutation.str.extract('(\d+)')
    one_mut_sortdf = onemut_df.join(pos_df, lsuffix='_caller', rsuffix='_other')\
        .set_index(0)\
        .fillna(0)\
        .sort_index()
    lowinput_df = one_mut_sortdf.loc[(one_mut_sortdf['input1']<=10) | (one_mut_sortdf['input2']<=10)] #Pandas dataframe filter with Multiple conditions
    one_mut_sortdf = one_mut_sortdf.loc[(one_mut_sortdf['input1']>10) & (one_mut_sortdf['input2']>10)]
    correlation = one_mut_sortdf.corr(method='pearson')
    correlation.to_csv('results/nep_mut_correlation.csv')

    Analyzed_onemut_df = cal_sum(one_mut_sortdf,'input','input1','input2')
    cal_enri(Analyzed_onemut_df,'output1','input1')
    cal_enri(Analyzed_onemut_df, 'output2', 'input2')
    cal_mean(Analyzed_onemut_df,'average_RF','output1_RF','output2_RF')
    cal_min(Analyzed_onemut_df,'min_RF','output1_RF','output2_RF')
    # Concatenate vertically
    Analyzed_onemut_df = pd.concat([Analyzed_onemut_df, lowinput_df], axis=0)
    Analyzed_onemut_df['mutation_type'] = 'missense'
    Analyzed_onemut_df['mutation_type'][Analyzed_onemut_df.NEP_Mutation.str.contains('silent')] = 'silent'
    Analyzed_onemut_df['mutation_type'][
            (~Analyzed_onemut_df.NEP_Mutation.str.contains('silent')) & Analyzed_onemut_df.NEP_Mutation.str.contains('_')] = 'nonsense'

    Analyzed_onemut_df.to_csv('results/Analyzed_onemut.csv')


if __name__ == "__main__":
  main()