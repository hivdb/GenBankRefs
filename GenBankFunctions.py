from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import re
import Levenshtein


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


def string_edit_dist(string1, string2):
    distance = Levenshtein.distance(string1, string2)
    #print("Levenshtein Distance:", distance)
    return distance

def process_author_field(text):
    text = text.replace(' and', ',')
    authors = text.split('.,')
    return tuple(authors)

# This can be converted into 1 or 2 lines
def compare_accession_lists(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    if set1 == set2:
        return 1
    else:
        return 0

# This is currently not being used because we are requiring the author lists
# to be identical to merge references
def compare_author_lists(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    if set1.issubset(set2) or set2.issubset(set1):
        return 1
    else:
        shared_set = set1 & set2
        pcnt_set1 = len(shared_set) / len(set1)
        pcnt_set2 = len(shared_set) / len(set1)
        max_pcnt = max(pcnt_set1, pcnt_set2)
        return max_pcnt


def process_accession_lists(df, logfile):
    file = open(logfile, 'w')
    accession_list = df['accession']
    close_lists = {}
    for i, item1 in enumerate(accession_list):
        close_matches = []
        for j, item2 in enumerate(accession_list):
            if i >= j:
                continue
            score = compare_accession_lists(item1, item2)
            if score == 1:
                if is_reference_genome(item1[0]):
                    continue
                close_matches.append(j)
                file.write(f'{i}: {item1}\n{j}:{item2}\n\n')
        if len(close_matches) >=1:
            close_lists[i] = close_matches

    list_of_shared_indexes = []
    for index, value in close_lists.items():
        shared_indexes = set()
        shared_indexes.add(index)
        shared_indexes.update(set(value))
        subset_flag = False
        for items in list_of_shared_indexes:
            if shared_indexes.issubset(items):
                subset_flag = True
        if subset_flag == True:
            continue
        list_of_shared_indexes.append(shared_indexes)
    print(f'''number with same accessions: {len(list_of_shared_indexes)}:   
          {list_of_shared_indexes}''')

    indexes_of_rows_to_be_dropped = [
        j
        for i in list_of_shared_indexes
        for j in i
    ]
    print('Dropped', len(indexes_of_rows_to_be_dropped))
    list_of_new_rows = []
    for item in list_of_shared_indexes:
        new_row = merge_refs_sharing_accessions(df, list(item))
        list_of_new_rows.append(new_row)

    df = df.drop(indexes_of_rows_to_be_dropped)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)

    return list_of_new_rows

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

    list_of_shared_indexes = []
    for index, value in close_lists.items():
        shared_indexes = set()
        shared_indexes.add(index)
        shared_indexes.update(set(value))
        subset_flag = False
        for items in list_of_shared_indexes:
            if shared_indexes.issubset(items):
                subset_flag = True
        if subset_flag == True:
            continue
        list_of_shared_indexes.append(shared_indexes)
    
    print(f'''number with shared author title: {len(list_of_shared_indexes)}:   
          {list_of_shared_indexes}''')        
    indexes_of_rows_to_be_dropped = [
        j
        for i in list_of_shared_indexes
        for j in i
    ]

    print('To be dropped', len(indexes_of_rows_to_be_dropped))
    list_of_new_rows = []
    for item in list_of_shared_indexes:
        new_row = merge_refs_sharing_author_title(df, list(item))
        list_of_new_rows.append(new_row)

    df = df.drop(indexes_of_rows_to_be_dropped)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)
    return list_of_new_rows


def merge_refs_sharing_accessions(df, shared_indexes):
    #print(df)
    new_row = {}
    author_list_list = df.loc[shared_indexes, 'author list'].tolist()
    authors_list = df.loc[shared_indexes, 'authors'].tolist()
    titles_list = df.loc[shared_indexes, 'title'].tolist()
    journal_list = df.loc[shared_indexes, 'journal'].tolist()
    pmid_list = df.loc[shared_indexes, 'pmid'].tolist()
    year_list = df.loc[shared_indexes, 'year'].tolist()
    print(year_list)
    accession_list = df.loc[shared_indexes, 'accession'].tolist()
    new_row['author list'] = [max(author_list_list, key=len)]
    new_row['authors'] = [max(authors_list, key=len)]
    new_row['year'] = [','.join(str(i) for i in year_list if i)]
    new_row['title'] = [max(titles_list, key=len)]
    new_row['pmid'] = ['; '.join(i for i in pmid_list if i)]
    new_row['journal'] = ['; '.join(journal_list)]
    new_row['accession'] = accession_list[:1]
    new_row_df = pd.DataFrame(new_row)
    return new_row_df

def merge_refs_sharing_author_title(df, shared_indexes):
    new_row = {}
    author_list_list = df.loc[shared_indexes, 'author list'].tolist()
    authors_list = df.loc[shared_indexes, 'authors'].tolist()
    titles_list = df.loc[shared_indexes, 'title'].tolist()
    journal_list = df.loc[shared_indexes, 'journal'].tolist()
    pmid_list = df.loc[shared_indexes, 'pmid'].tolist()
    year_list = df.loc[shared_indexes, 'year'].tolist()
    accession_list = df.loc[shared_indexes, 'accession'].tolist()
    new_row['author list'] = author_list_list[:1] 
    new_row['authors'] = authors_list[:1]
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

# def fetch_popset_info(uid):
#     # Query PopSet database with the UID
#     handle = Entrez.efetch(db="popset", id=uid, rettype="acc", retmode="text")
#     record = Entrez.read(handle)
#     print(record)
#     acc_ids = [id.strip() for id in handle]
#     handle.close()

# def search_popsets_for_virus(virus_name):
#     # Perform a search in the PopSet database using the virus name
#     search_query = f"{virus_name}[Organism]"
#     handle = Entrez.esearch(db="popset", term=search_query, retmax=200)
#     record = Entrez.read(handle)
#     handle.close()

#     # Get the list of UIDs (PopSet IDs) from the search results
#     uid_list = record["IdList"]
#     print(uid_list)

#     # If there are results, fetch detailed info about each PopSet
#     if uid_list:
#         print(f"Found {len(uid_list)} PopSets for {virus_name}.")
#         for uid in uid_list:
#             fetch_popset_info(uid)  # Function from previous example
#     else:
#         print(f"No PopSets found for {virus_name}.")

