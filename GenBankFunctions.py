from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import re
import Levenshtein


Entrez.email = "rshafer.stanford.edu"

def extract_year_from_journal(text):
    match = re.search(r'\((\d{4})\)', text)
    if match:
        #print(match.group(1))
        return match.group(1)
    match = re.search (r'\d{2}-[A-Z]{3}-(\d{4})', text)
    if match:
        #print(match.group(1))
        return match.group(1)
    return ''

def string_edit_dist(string1, string2):
    distance = Levenshtein.distance(string1, string2)
    print("Levenshtein Distance:", distance)
    return distance

def process_author_field(text):
    text = text.replace(' and', ',')
    authors = text.split('.,')
    return tuple(authors)

def compare_accession_lists(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    if set1 == set2:
        return 1
    else:
        return 0

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



def process_accession_lists(df):
    new_df_consolidate_rows_with_same_accessions = df
    accession_list = df['accession']
    close_lists = {}
    for i, item1 in enumerate(accession_list):
        close_matches = []
        for j, item2 in enumerate(accession_list):
            if i >= j:
                continue
            score = compare_accession_lists(item1, item2)
            if score == 1:
                if item1[0].startswith('NM') or item1[0].startswith('NC') or item1[0].startswith('NG'):
                    continue
                close_matches.append(j)
                #print("Accesion list: ", item1, "\nAccession list: ", item2, "\n", score, "\n")
        if len(close_matches) >=1:
            close_lists[i] = close_matches

    list_of_indexes_with_same_accessions = []
    for index, value in close_lists.items():
        indexes_with_same_accessions = set()
        indexes_with_same_accessions.add(index)
        indexes_with_same_accessions.update(set(value))
        subset_flag = False
        for items in list_of_indexes_with_same_accessions:
            if indexes_with_same_accessions.issubset(items):
                subset_flag = True
        if subset_flag == True:
            continue
        list_of_indexes_with_same_accessions.append(indexes_with_same_accessions)
    print('same accessions', len(list_of_indexes_with_same_accessions))

    indexes_of_rows_to_be_dropped = [
        j
        for i in list_of_indexes_with_same_accessions
        for j in i
    ]
    print('Dropped', len(indexes_of_rows_to_be_dropped))
    list_of_new_rows = []
    for item in list_of_indexes_with_same_accessions:
        new_row = merge_refs_sharing_accessions(df, list(item))
        list_of_new_rows.append(new_row)

    df = df.drop(indexes_of_rows_to_be_dropped)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)

    return list_of_new_rows



def process_author_sets(author_list_column):
    close_lists = {}
    for i, item1 in enumerate(author_list_column):
        for j, item2 in enumerate(author_list_column):
        # Skip comparison if the indices are the same
            if i >= j:
                continue
            score = compare_author_lists(item1, item2)
            if (score < 0.8):
                continue
            close_lists[(i,j)] = score
            print("Authorlist1: ", item1, "\nAuthorList2: ", item2, "\n", score, "\n")

def merge_refs_sharing_accessions(df, indexes_with_same_accessions):
    new_row = {}
    author_list_list = df.loc[indexes_with_same_accessions, 'author list'].tolist()
    authors_list = df.loc[indexes_with_same_accessions, 'authors'].tolist()
    titles_list = df.loc[indexes_with_same_accessions, 'title'].tolist()
    journal_list = df.loc[indexes_with_same_accessions, 'journal'].tolist()
    pmid_list = df.loc[indexes_with_same_accessions, 'pmid'].tolist()
    year_list = df.loc[indexes_with_same_accessions, 'year'].tolist()
    accession_list = df.loc[indexes_with_same_accessions, 'accession'].tolist()
    new_row['author list'] = [max(author_list_list, key=len)]
    new_row['authors'] = [max(authors_list, key=len)]
    new_row['year'] = [','.join(i for i in year_list if i)]
    new_row['title'] = [max(titles_list, key=len)]
    new_row['pmid'] = ['; '.join(i for i in pmid_list if i)]
    new_row['journal'] = ['; '.join(journal_list)]
    new_row['accession'] = accession_list[:1]
    new_row_df = pd.DataFrame(new_row)
    # print(new_row_df, "\n\n")

    return new_row_df




def fetch_popset_info(uid):
    # Query PopSet database with the UID
    handle = Entrez.efetch(db="popset", id=uid, rettype="acc", retmode="text")
    record = Entrez.read(handle)
    print(record)
    acc_ids = [id.strip() for id in handle]
    handle.close()

def search_popsets_for_virus(virus_name):
    # Perform a search in the PopSet database using the virus name
    search_query = f"{virus_name}[Organism]"
    handle = Entrez.esearch(db="popset", term=search_query, retmax=200)
    record = Entrez.read(handle)
    handle.close()

    # Get the list of UIDs (PopSet IDs) from the search results
    uid_list = record["IdList"]
    print(uid_list)

    # If there are results, fetch detailed info about each PopSet
    if uid_list:
        print(f"Found {len(uid_list)} PopSets for {virus_name}.")
        for uid in uid_list:
            fetch_popset_info(uid)  # Function from previous example
    else:
        print(f"No PopSets found for {virus_name}.")

