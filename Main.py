from Bio import SeqIO
from Bio import Entrez
import pandas as pd


Entrez.email = "rshafer.stanford.edu"

pd.set_option('display.max_rows', 100)

genbank_file = "CCHF.gb"

from GenBankFunctions import extract_year_from_journal, process_author_field, process_accession_lists, process_authors_titles

def extract_organism(annotations):
    organism = annotations.get("source")
    #print("\nOrganism:", organism, "\n")

def extract_references(annotations, accession):
    ref_list = annotations.get("references")
    ref_data = []
    for ref in ref_list:
        #print(dir(ref))
        ref_items = {}
        #print ("\nReference:", ref, "\n")
        ref_items["accession"] = accession
        ref_items["authors"] = ref.authors
        ref_items["title"] = ref.title
        ref_items["journal"] = ref.journal
        ref_items["pmid"] = ref.pubmed_id
        ref_data.append(ref_items)
    return ref_data

def extract_features(features, accession):
    #feature_data = []
    feature_items = {}
    feature_items["accession"] = accession
    for feature in features:
        #print(f"\nFeature type: {feature.type}\n")
        for key, value in feature.qualifiers.items():
            if key == "translation":
                continue
            key = key + '_' + feature.type
            feature_items[key] = value[0]
            #print(f"Key {key}: value {value}\n")
    #feature_data.append(feature_items)
    return feature_items

# Example usage
#search_popsets_for_virus("Orthonairovirus")

reference_list = []
feature_list = []
with open(genbank_file, "r") as handle:
    count=0
    for record in SeqIO.parse(handle, "genbank"):
        count +=1
        #print(f"ID: {record.id}")
        #print(f"Description: {record.description}")
        #print(f"Sequence Length: {len(record.seq)}")
        #print(f"Annotations: {record.annotations}")
        #print(f"Features (length): {len(record.features)}")
        #print(f"Features: {record.features}")
        #print("-" * 40)  # Separator between entries
        #extract_organism(record.annotations)
        ref_data = extract_references(record.annotations, record.id)
        reference_list.extend(ref_data)
        feature_data = extract_features(record.features, record.id)
        feature_list.append(feature_data)
        #print(feature_data)
        #extract_features(record.features)
        #if count > 5:
        #    break

## Aggregate by reference
reference_df = pd.DataFrame(reference_list)
reference_df['year'] = reference_df['journal'].apply(extract_year_from_journal)
reference_df['year'] = pd.to_numeric(reference_df['year'], errors='coerce')
reference_df['journal'] = reference_df['journal'].str.replace(r"Submitted \(\d{2}-[A-Z]{3}-\d{4}\)", "", regex=True)
reference_df['journal'] = reference_df['journal'].str.replace(r"(Patent).*", r"\1", regex=True)
reference_df['author list'] = reference_df['authors'].apply(process_author_field)

print(len(reference_df))
grouped_ref_df = reference_df.groupby(['authors', 'author list', 'title', 'journal', 'pmid', 'year'])['accession'].apply(list).reset_index()
print(len(grouped_ref_df))
#print(grouped_ref_df.head(100))
grouped_ref_df.to_excel("Grouped_Refs.xlsx")

merged_ref_df = process_accession_lists(grouped_ref_df, "accession_code.txt")
print(len(merged_ref_df))
merged_ref_df.to_excel("Merged_Accessions.xlsx")

merged_ref_df = process_authors_titles(merged_ref_df, "authors_title_code.txt") 
merged_ref_df.to_excel("Merged_Author_Titles.xlsx")
#process_author_sets(grouped_ref_df['author list'])


#reference_df.to_csv("CCHF_GenBankRefs.csv", index=False)
#df_refs = pd.read_csv('CCHF_GenBankRefs.csv')
# Create a new dataframe containing just the 2nd to 5th columns
# new_ref_df = reference_df.iloc[:, 1:5]
# unique_refs_df = new_ref_df.drop_duplicates()
# print(unique_refs_df)
# print(len(unique_refs_df))


# Add a reference back to the original "accession" column from the original dataframe
# We will do this by reattaching the first column ('accession') to the unique dataframe
#unique_df_with_accession = pd.concat([df['accession'], new_df], axis=1).drop_duplicates()

#import ace_tools as tools; tools.display_dataframe_to_user(name="Unique DataFrame with Accession", dataframe=unique_df_with_accession)






# features_df = pd.DataFrame(feature_list)
#features_df.to_csv("CCHF_GenBankFeatures.csv", index = False)

# value_counts = features_df["organism_source"].value_counts()
# print(value_counts)
# value_counts = features_df["mol_type_source"].value_counts()
# print(value_counts)
# value_counts = features_df["strain_source"].value_counts()
# print(value_counts)
# value_counts = features_df["db_xref_source"].value_counts()
# print(value_counts)
# value_counts = features_df["segment_source"].value_counts()
# print(value_counts)
# value_counts = features_df["geo_loc_name_source"].value_counts()
# print(value_counts)
# value_counts = features_df["product_CDS"].value_counts()
# print(value_counts)










