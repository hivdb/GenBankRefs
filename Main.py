from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import ast

Entrez.email = "rshafer.stanford.edu"
pd.set_option('display.max_rows', 100)
genbank_file = "CCHF.gb"

from GenBankFunctions import extract_year_from_journal, process_author_field, process_accession_lists, process_authors_titles

#def extract_organism(annotations):
#    organism = annotations.get("source")
#    #print("\nOrganism:", organism, "\n")

def extract_references(annotations, accession):
    ref_list = annotations.get("references")
    ref_data = []
    for ref in ref_list:
        ref_items = {}
        ref_items["accession"] = accession
        ref_items["authors"] = ref.authors
        ref_items["title"] = ref.title
        ref_items["journal"] = ref.journal
        ref_items["pmid"] = ref.pubmed_id
        ref_data.append(ref_items)
    return ref_data

#columns_of_interest = []
def extract_features(features, accession):
    feature_items = {}
    feature_items["accession"] = accession
    for feature in features:
        for key, value in feature.qualifiers.items():
            if key == "translation":
                continue
            key = key + '_' + feature.type
            feature_items[key] = value[0]
           #print(f"Key {key}: value {value}\n")
    return feature_items


reference_list = []
feature_list = []
with open(genbank_file, "r") as handle:
    count=0
    for record in SeqIO.parse(handle, "genbank"):
        count +=1
        acc_num = record.id
        seq_len = len(record.seq)
        description = record.description
        record_date = record.annotations['date']
        organism = record.annotations['organism']
     
        ref_data = extract_references(record.annotations, record.id)
        reference_list.extend(ref_data)
        
        feature_data = extract_features(record.features, record.id)
        #print("\nFeature data:", feature_data, "\n")

        segment_source = feature_data.get('segment_source', '')
        cds = feature_data.get('product_CDS', '')
        host = feature_data.get('host_source', '')
        isolate_source = feature_data.get('isolation_source_source', '')
        isolate_source_name = feature_data.get('isolate_source', '')
        country_region = feature_data.get('geo_loc_name_source', '')
        collection_date = feature_data.get('collection_date_source', '')

        print(f"Acc_num: {acc_num}")
        print(f"Organism: {organism}")
        print(f"Description: {description}")
        print(f"Record date: {record_date}")
        print(f"Sequence Length: {seq_len}")
        print("Segment source:", segment_source)
        print("CDS:", cds)
        print("Host:", host)
        print("Isolation source:", isolate_source)
        print("Isolate_name:", isolate_source_name)
        print("Country region:", country_region)
        print("Collection date:", collection_date)

        feature_list.append(feature_data)
        print("___________________________________________________")
        if count > 10:
            break

## Aggregate by reference
excluded_accessions = ['NM_010185.4', 'NM_010508.2', 'NM_134350.2', 'NM_021268.2', 'NM_009283.4', \
                       'NM_001205313.1', 'NM_001205314.1', 'NM_001357627.1', 'NM_009283.4', \
                       'NR_104124.3', 'NR_104125.3', 'NR_186227.2', 'NR_186228.2', 'NR_186229.2', \
                       'NR_186230.2', 'NR_186231.2', 'NR_126359.1', \
                       'JAHWGI010000070.1', 'JAHWGI010000979.1', 'JAHWGI010001134.1', 'JAHWGI010001142.1', 'JAHWGI010001356.1']
reference_df = pd.DataFrame(reference_list)
reference_df['year'] = reference_df['journal'].apply(extract_year_from_journal)
reference_df['year'] = pd.to_numeric(reference_df['year'], errors='coerce')
reference_df['year'] = reference_df['year'].apply(lambda x: '' if pd.isna(x) else int(x))
reference_df['journal'] = reference_df['journal'].str.replace(r"Submitted \(\d{2}-[A-Z]{3}-\d{4}\)", "", regex=True)
reference_df['journal'] = reference_df['journal'].str.replace(r"(Patent).*", r"\1", regex=True)
reference_df['authors'] = reference_df['authors'].apply(process_author_field)
reference_df = reference_df[~reference_df['accession'].isin (excluded_accessions)]
print("Number of original entries: ", len(reference_df))

grouped_ref_df = reference_df.groupby(['authors', 'title', 'journal', 'pmid', 'year'])['accession'].apply(list).reset_index()
grouped_ref_df['accession'] = grouped_ref_df['accession'].apply(lambda x: ', '.join(x))
print("Number of entries following aggregation by metadata: ", len(grouped_ref_df))
grouped_ref_df.to_excel("CCHF_Grouped_Refs.xlsx")

merged_ref_df = process_accession_lists(grouped_ref_df)
print("Number of entries following aggregation by accession numbers: ", len(merged_ref_df))
merged_ref_df.to_excel("CCHF_Merged_Accessions.xlsx")

merged_ref_df = process_authors_titles(merged_ref_df) 
print("Number of entries following aggregation by metadata: ", len(merged_ref_df))
merged_ref_df.to_excel("CCHF_Merged_Author_Titles.xlsx")

features_df = pd.DataFrame(feature_list)
features_df.to_csv("CCHF_GenBankFeatures.csv", index = False)

# Add a reference back to the original "accession" column from the original dataframe
# We will do this by reattaching the first column ('accession') to the unique dataframe
#unique_df_with_accession = pd.concat([df['accession'], new_df], axis=1).drop_duplicates()
#import ace_tools as tools; tools.display_dataframe_to_user(name="Unique DataFrame with Accession", dataframe=unique_df_with_accession)













