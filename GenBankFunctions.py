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
    accession_list = df('accession')
    close_lists = {}
    for i, item1 in enumerate(accession_list):
        for j, item2 in enumerate(accession_list):
        # Skip comparison if the indices are the same
            score = compare_accession_lists(item1, item2)
            if score == 1:
                close_lists[(i,j)] = score
            print("Authorlist1: ", item1, "\nAuthorList2: ", item2, "\n", score, "\n")


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

