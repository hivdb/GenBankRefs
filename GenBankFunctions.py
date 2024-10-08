from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import re


Entrez.email = "rshafer.stanford.edu"

def extract_year_from_journal(text):
    match = re.search(r'\((\d{4})\)', text)
    if match:
        return match.group(1)
    match = re.search (r'\d{2}-[A-Z]{3}-(\d{4})', text)
    if match:
        return match.group(1)
    return ''

def is_reference_genome(acc):
    prefix_list = ['NC', 'NG','NM', 'NR']
    for item in prefix_list:
        if acc.startswith(item):
            return True
    return False


def process_author_field(text):
    text = text.replace(' and', ',')
    authors = text.split('.,')
    return tuple(authors)


# The score is 1 if the author list is the same and 1 of the following criteria is met
# (1) One or both of the titles are 'Direct Submission" and the years are within 1. If the years
# are not present they are treated as being the same
# (2) The titles are the same (and not both 'Direct Submission')
def compare_authors_titles(row_i, row_j):
    score = 0
    author_list_i = row_i['author list']
    author_list_j = row_j['author list']
    title_i = row_i['title']
    title_j = row_j['title']
    year_i = row_i['year']
    year_j = row_j['year']
  
    if author_list_i == author_list_j and len(author_list_i) > 1:       
        if title_i == 'Direct Submssion' or title_j == 'Direct Submission':
            if pd.notna(year_i) and pd.notna(year_j):
                if year_i == year_j:
                    score = 1
                else:
                    score = 0
            else:
                score = 1     
        elif title_i == title_j:
            score = 1             
    return score


#Convert {key: [values], key: [values]] to a list of sets comprising {key, values}
#In this program, each key is the index of a row and the values contain
# one or more indexes of rows that share some property
#The subset flag code ensures sets do not share items with one another
def convert_dict_to_list_of_sets(dict):
    list_of_sets = []
    list_of_items = []
    for key, values in dict.items():
        key_plus_values = set()
        key_plus_values.add(key)
        key_plus_values.update(set(values))    
        subset_flag = False
        for items in list_of_sets:
            if key_plus_values.issubset(items):
                subset_flag = True
        if subset_flag == True:
            continue
        list_of_sets.append(key_plus_values)
    
    list_of_items = [
        j
        for i in list_of_sets
        for j in i
    ]

    return (list_of_sets, list_of_items)


def process_accession_lists(df, logfile):
    file = open(logfile, 'w')
    accession_list = df['accession']
    close_lists = {}
    for i, item1 in enumerate(accession_list):
        close_matches = []
        for j, item2 in enumerate(accession_list):
            if i >= j:
                continue
            if set(item1) == set(item2) and is_reference_genome(item1[0]) == False:
                close_matches.append(j)
                file.write(f'{i}: {item1}\n{j}:{item2}\n\n')
        if len(close_matches) >=1:
            close_lists[i] = close_matches

    list_of_sets_w_shared_indexes, complete_list_of_shared_indexes = convert_dict_to_list_of_sets(close_lists)  
    print(f'''No with shared accessions: {len(list_of_sets_w_shared_indexes)}: {list_of_sets_w_shared_indexes}''')
    print(f'''To be dropped: {len(complete_list_of_shared_indexes)}: {complete_list_of_shared_indexes}''')
    
    list_of_new_rows = []
    for item in list_of_sets_w_shared_indexes:
        new_row = merge_refs_sharing_accessions(df, list(item))
        list_of_new_rows.append(new_row)

    df = df.drop(complete_list_of_shared_indexes)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)
    return list_of_new_rows


def process_authors_titles(df, logfile):
    file = open(logfile, 'w')
    close_lists = {} 
    for i, row_i in df.iterrows():
        close_matches = []
        for j, row_j in df.iterrows(): 
            if i >= j:
                continue
            score = compare_authors_titles(row_i, row_j)
            if score == 1:
                close_matches.append(j)
                file.write(f'{i}: {row_i}\n {j}: {row_j}\n\n')
            if len(close_matches) >=1:
                close_lists[i] = close_matches
    print("Close lists:", close_lists)

    list_of_sets_w_shared_indexes, complete_list_of_shared_indexes = convert_dict_to_list_of_sets(close_lists)  
    print(f'''No with shared accessions: {len(list_of_sets_w_shared_indexes)}: {list_of_sets_w_shared_indexes}''')
    print(f'''To be dropped: {len(complete_list_of_shared_indexes)}: {complete_list_of_shared_indexes}''')      

    list_of_new_rows = []
    for item in list_of_sets_w_shared_indexes:
        new_row = merge_refs_sharing_author_title(df, list(item))
        list_of_new_rows.append(new_row)

    df = df.drop(complete_list_of_shared_indexes)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)
    return list_of_new_rows


def merge_refs_sharing_accessions(df, shared_indexes):
    new_row = {}
    author_list = df.loc[shared_indexes, 'author list'].tolist()
    authors = df.loc[shared_indexes, 'authors'].tolist()
    titles_list = df.loc[shared_indexes, 'title'].tolist()
    journal_list = df.loc[shared_indexes, 'journal'].tolist()
    pmid_list = df.loc[shared_indexes, 'pmid'].tolist()
    year_list = df.loc[shared_indexes, 'year'].tolist()
    accession_list = df.loc[shared_indexes, 'accession'].tolist()
    new_row['author list'] = [max(author_list, key=len)]
    new_row['authors'] = [max(authors, key=len)]
    new_row['year'] = [','.join(str(i) for i in year_list if i)]
    new_row['title'] = [max(titles_list, key=len)]
    new_row['pmid'] = ['; '.join(i for i in pmid_list if i)]
    new_row['journal'] = ['; '.join(journal_list)]
    new_row['accession'] = accession_list[:1]
    new_row_df = pd.DataFrame(new_row)
    return new_row_df

def merge_refs_sharing_author_title(df, shared_indexes):
    new_row = {}
    author_list = df.loc[shared_indexes, 'author list'].tolist()
    authors = df.loc[shared_indexes, 'authors'].tolist()
    titles_list = df.loc[shared_indexes, 'title'].tolist()
    journal_list = df.loc[shared_indexes, 'journal'].tolist()
    pmid_list = df.loc[shared_indexes, 'pmid'].tolist()
    year_list = df.loc[shared_indexes, 'year'].tolist()
    accession_list = df.loc[shared_indexes, 'accession'].tolist()
    new_row['author list'] = author_list[:1] 
    new_row['authors'] = authors[:1]
    new_row['year'] = [';'.join(str(i) for i in year_list if i)]
    new_row['title'] = ['\n'.join(i for i in titles_list)]
    new_row['pmid'] = ['; '.join(i for i in pmid_list if i)]
    new_row['journal'] = ['; '.join(journal_list)]
    combined_accessions = set()
    for accessions in accession_list:
        combined_accessions.update(set(accessions))
    combined_list_accessions = list(combined_accessions)
    combined_list_accessions.sort()
    new_row['accession'] = ','.join(combined_list_accessions)
    new_row_df = pd.DataFrame(new_row)
    return new_row_df

