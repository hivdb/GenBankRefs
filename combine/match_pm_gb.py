import pandas as pd
from collections import defaultdict
import Levenshtein


def match_pm_gb(pubmed, genbank):

    match_by_title = []
    match_by_pmid = []
    match_by_acc = []
    matched_pubmed_indices = []
    genbank_unmatch = []

    for index, row in genbank.iterrows():
        pmid = row['pmid']

        result = pubmed[pubmed['PMID'] == pmid]

        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            match_by_pmid.append([row, result])
            continue

        title = row['title']

        result = pubmed[pubmed['Title'].apply(
            lambda x: Levenshtein.distance(x.lower(), title.lower()) < 5)]

        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            match_by_title.append([row, result])
            continue

        accession_list = row['accession']
        accession_prefix_list = set([
            a.strip()[:6]
            for a in accession_list.split(',')
        ])

        # if 'JX908640' in accession_list:
        #     print(accession_prefix_list)

        result = search_access_prefix(pubmed, accession_prefix_list)
        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            match_by_acc.append([row, result])
            continue

        genbank_unmatch.append(row)

    genbank_match = match_by_title + match_by_pmid + match_by_acc

    print("Genbank match by pmid", len(match_by_pmid))
    print("Genbank match by title", len(match_by_title))
    print("Genbank match by acc", len(match_by_acc))
    print('Genbank Matched:', len(genbank_match))

    pubmed_match = match_pubmed2genbank(genbank_match)

    print("Pubmed match by pmid", len(match_pubmed2genbank(match_by_pmid)))
    print("Pubmed match by title", len(match_pubmed2genbank(match_by_title)))
    print("Pubmed match by acc", len(match_pubmed2genbank(match_by_acc)))
    print('Pubmed Matched:', len(pubmed_match))

    pubmed_unmatch = pubmed.drop(index=matched_pubmed_indices)

    return pubmed_match, pubmed_unmatch, genbank_unmatch


def search_access_prefix(pubmed, accession_prefix_list):
    found = pd.DataFrame()

    for acc_prefix in accession_prefix_list:
        result = pubmed[pubmed['GenBank'].apply(
            lambda x: acc_prefix.lower() in str(x).lower())]
        if not result.empty:
            found = pd.concat([found, result])

    return found


def match_pubmed2genbank(genbank_match):
    pubmed_match = defaultdict(dict)
    for g, pubmed_list in genbank_match:
        for r, i in pubmed_list.iterrows():
            pubmed_match[r]['pubmed'] = i
            if 'genbank' not in pubmed_match[r]:
                pubmed_match[r]['genbank'] = []

            pubmed_match[r]['genbank'].append(g)

    pubmed_match = [
        [v['pubmed'], v['genbank']]
        for v in pubmed_match.values()
    ]

    return pubmed_match
