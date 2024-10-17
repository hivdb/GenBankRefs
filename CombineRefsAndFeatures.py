import pandas as pd
import re

from GenBankFunctions import (merge_feature_rows)

pd.set_option('display.max_rows', 100)

ref_df = pd.read_excel('CCHF_Merged_Author_Titles_Oct15.xlsx')
features_df = pd.read_excel('CCHF_GenBankFeatures_Oct15.xlsx')

#print(ref_df)
#print(features_df)

combined_df = ref_df.copy()

feature_columns = ['Organisms', 'RecordYears',  'Hosts', 'Countries', 
                   'IsolateYears', 'CDS', 'SeqLens', 'AlignLens', 'PcntIDs']

combined_df[feature_columns] = 'None'
print(combined_df)

count = 0
for index, row in combined_df.iterrows():
    count += 1
    accession_string = row['accession']
    accession_list = accession_string.split(', ')
    print("\n", index)
    print(accession_list)
    features_rows = features_df[features_df['acc_num'].isin(accession_list)]
    num_rows = len(features_rows)
    new_dict = merge_feature_rows(features_rows)
    print(new_dict)
    combined_df.at[index, 'Organisms'] = new_dict['Organisms']
    combined_df.at[index, 'RecordYears'] = new_dict['RecordYears']
    combined_df.at[index, 'Hosts'] = new_dict['Hosts']
    combined_df.at[index, 'Countries'] = new_dict['Countries']
    combined_df.at[index, 'IsolateYears'] = new_dict['IsolateYears']
    combined_df.at[index, 'CDS'] = new_dict['CDS']
    combined_df.at[index, 'SeqLens'] = new_dict['SeqLens']
    combined_df.at[index, 'AlignLens'] = new_dict['AlignLens']
    combined_df.at[index, 'PcntIDs'] = new_dict['PcntIDs']

combined_df.to_csv("CCHF_GenBank_Combined.csv", index = False) 

