import pandas as pd
from collections import defaultdict
import Levenshtein


def match_pm_gb(pubmed, genbank):

    match_by_title_list = []
    match_by_pmid_list = []
    match_by_acc_list = []
    matched_pubmed_indices = []
    genbank_unmatch_list = []

    for index, row in genbank.iterrows():
        pmid = row['PMID']

        result = pubmed[pubmed['PMID'] == pmid]

        match_by_acc = None
        match_by_pmid = None
        match_by_title = None

        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            match_by_pmid = [row, result, index]
            match_by_pmid_list.append(match_by_pmid)
            # continue

        title = row['Title'].replace('Direct Submission,', '').replace(', Direct Submission', '').strip()

        # Pubmed title always exists
        result = pubmed[pubmed['Title'].apply(
            lambda x: Levenshtein.distance(x.lower(), title.lower()) < 5)]

        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            match_by_title = [row, result, index]
            match_by_title_list.append(match_by_title)

        accession_list = row['accession']
        accession_prefix_list = set([
            a.strip()[:6]
            for a in accession_list.split(',')
        ])

        result = search_access_prefix(pubmed, accession_prefix_list)
        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            match_by_acc = [row, result, index]
            match_by_acc_list.append(match_by_acc)

        if match_by_acc or match_by_pmid or match_by_title:
            pass
        else:
            genbank_unmatch_list.append(row)

    genbank_match_list = match_by_title_list + match_by_pmid_list + match_by_acc_list

    print("Genbank match by pmid", len(set(i[-1] for i in match_by_pmid_list)))
    print("Genbank match by title", len(set(i[-1] for i in match_by_title_list)))
    print("Genbank match by acc", len(set(i[-1] for i in match_by_acc_list)))
    print('Genbank Matched:', len(set(i[-1] for i in genbank_match_list)))

    pubmed_match = match_pubmed2genbank(genbank_match_list)

    genbank_match = {}
    for row, result, index in genbank_match_list:
        genbank_match[index] = row

    genbank_match = genbank_match.values()

    print("Pubmed match by pmid", len(match_pubmed2genbank(match_by_pmid_list)))
    print("Pubmed match by title", len(match_pubmed2genbank(match_by_title_list)))
    print("Pubmed match by acc", len(match_pubmed2genbank(match_by_acc_list)))
    print('Pubmed Matched:', len(pubmed_match))

    pubmed_unmatch = pubmed.drop(index=matched_pubmed_indices)

    return pubmed_match, pd.DataFrame(genbank_match), pubmed_unmatch, pd.DataFrame(genbank_unmatch_list)


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
    for g, pubmed_list, index in genbank_match:
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
