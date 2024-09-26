from Bio import SeqIO
import pandas as pd

genbank_file = "CCHF.gb"


def extract_organism(annotations):
    organism = annotations.get("source")
    #print("\nOrganism:", organism, "\n")

def extract_references(annotations, accession):
    ref_list = annotations.get("references")
    ref_data = []
    for ref in ref_list:
        #print(dir(ref))
        ref_items = {}
        print ("\nReference:", ref, "\n")
        ref_items["accession"] = accession
        ref_items["authors"] = ref.authors
        ref_items["title"] = ref.title 
        ref_items["journal"] = ref.journal
        ref_items["pmid"] = ref.pubmed_id
        ref_data.append(ref_items)
    return ref_data

def extract_features(features):
    for feature in features:
        #print(f"\nFeature type: {feature.type}\n")
        for key, value in feature.qualifiers.items():
            print(f"Key {key}: value {value}\n")

reference_list = []
with open(genbank_file, "r") as handle:
    count=0
    for record in SeqIO.parse(handle, "genbank"):
        count +=1
        print(f"ID: {record.id}")
        print(f"Description: {record.description}")
        print(f"Sequence Length: {len(record.seq)}")
        #print(f"Annotations: {record.annotations}")
        #print(f"Features (length): {len(record.features)}")
        #print(f"Features: {record.features}")
        #print("-" * 40)  # Separator between entries
        #extract_organism(record.annotations)
        ref_data = extract_references(record.annotations, record.id)
        #extract_features(record.features)
        reference_list.extend(ref_data)
        #if count > 5:
        #    break

reference_df = pd.DataFrame(reference_list)
reference_df.to_csv("CCHF_GenBankRefs.csv", index=False)







