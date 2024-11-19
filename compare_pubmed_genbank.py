from openai import OpenAI
import pandas as pd
from dataclasses import dataclass, asdict
from pathlib import Path
from tqdm import tqdm
from time import sleep
from collections import defaultdict


def search_access_prefix(pubmed, accession_prefix_list):
    found = pd.DataFrame()

    for acc_prefix in accession_prefix_list:
        result = pubmed[pubmed['GenBank'].apply(lambda x: acc_prefix.lower() in str(x).lower())]
        if not result.empty:
            found = pd.concat([found, result])

    return found


def main():
    genbank_file = 'CCHF_Combined_10_29_saved.xlsx'
    pubmed_file = 'ReferenceSummary_Nov16.xlsx'

    genbank = pd.read_excel(genbank_file)
    pubmed = pd.read_excel(pubmed_file)

    genbank_match = []
    matched_pubmed_indices = []
    genbank_unmatch = []

    for index, row in genbank.iterrows():
        pmid = row['pmid']

        result = pubmed[pubmed['PMID'] == pmid]

        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            genbank_match.append([row, result])
            continue

        title = row['title']

        result = pubmed[pubmed['Title'].apply(lambda x: x.lower() in title.lower())]

        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            genbank_match.append([row, result])
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
            genbank_match.append([row, result])
            continue

        genbank_unmatch.append(row)

    pubmed.to_excel('Pubmed.xlsx')

    print('Genbank Matched:', len(genbank_match))

    pubmed_match = match_pubmed2genbank(genbank_match)
    print('Pubmed Matched:', len(pubmed_match))

    pubmed_unmatch = pubmed.drop(index=matched_pubmed_indices)

    # pubmed_unmatch.to_excel('pubmed_unmatch.xlsx')
    # pd.DataFrame(genbank_unmatch).to_excel('genbank_unmatch.xlsx')

    combined, columns = combine(pubmed_match, pubmed_unmatch, genbank_unmatch)
    combined.to_excel('combined.xlsx', index=False, columns=columns)


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


def combine(pubmed_match, pubmed_unmatch, genbank_unmatch):

    result = []
    for pubmed, genbank_list in pubmed_match:
        row = {
            'Authors': pubmed['Authors'],
            'Title': pubmed['Title'],
            'Journal': pubmed['Journal'],
            'Year': pubmed['Year'],
            'PMID': pubmed['PMID'],
            'Viruses (AI)': pubmed['Viruses'],
            'NumSeqs (AI)': pubmed['NumSeqs'],
            'Hosts (AI)': pubmed['Host'],
            'Specimen (AI)': pubmed['IsolateType'],
            'SampleYr (AI)': pubmed['SampleYr'],
            'Countries (AI)': pubmed['Country'],
            # 'Genes (AI)':
            'SeqMethod (AI)': pubmed['SeqMethod'],
            'CloneMethod (AI)': pubmed['CloneMethod'],
            'GenBank (AI)': pubmed['GenBank'],

            'Viruses (GB)': ';'.join([genbank['Organisms'] for genbank in genbank_list]),
            'NumSeqs (GB)': len([
                i
                for genbank in genbank_list
                for i in genbank['accession'].split(',')
            ]),
            'Hosts (GB)': ';'.join([genbank['Hosts'] for genbank in genbank_list]),
            'Specimen (GB)': '',
            'SampleYr (GB)':  ';'.join([genbank['IsolateYears'] for genbank in genbank_list]),
            'Countries (GB)': ';'.join([genbank['Countries'] for genbank in genbank_list]),
            # 'Genes (GB)':
            'SeqMethod (GB)': '',
            'CloneMethod (GB)': '',
            'GenBank (GB)': ';'.join([genbank['accession'] for genbank in genbank_list]),
            'AlignLens (GB)': ';'.join([genbank['AlignLens'] for genbank in genbank_list]),
            'PcntIDs (GB)': ';'.join([str(genbank['PcntIDs']) for genbank in genbank_list]),
        }

        result.append(row)

    for r, pubmed in pubmed_unmatch.iterrows():
        row = {
            'Authors': pubmed['Authors'],
            'Title': pubmed['Title'],
            'Journal': pubmed['Journal'],
            'Year': pubmed['Year'],
            'PMID': pubmed['PMID'],
            'Viruses (AI)': pubmed['Viruses'],
            'NumSeqs (AI)': pubmed['NumSeqs'],
            'Hosts (AI)': pubmed['Host'],
            'Specimen (AI)': pubmed['IsolateType'],
            'SampleYr (AI)': pubmed['SampleYr'],
            'Countries (AI)': pubmed['Country'],
            # 'Genes (AI)':
            'SeqMethod (AI)': pubmed['SeqMethod'],
            'CloneMethod (AI)': pubmed['CloneMethod'],
            'GenBank (AI)': pubmed['GenBank'],
        }

        result.append(row)

    for genbank in genbank_unmatch:
        numSeqs = len(genbank['accession'].split(','))
        row = {
            'Authors': genbank['authors'],
            'Title': genbank['title'],
            'Journal': genbank['journal'],
            'Year': genbank['year'],
            'PMID': genbank['pmid'],
            'Viruses (GB)': genbank['Organisms'],
            'NumSeqs (GB)': numSeqs,
            'Hosts (GB)': genbank['Hosts'],
            'Specimen (GB)': '',
            'SampleYr (GB)': genbank['IsolateYears'],
            'Countries (GB)': genbank['Countries'],
            # 'Genes (GB)':
            'SeqMethod (GB)': '',
            'CloneMethod (GB)': '',
            'GenBank (GB)': genbank['accession'],
            'AlignLens (GB)': genbank['AlignLens'],
            'PcntIDs (GB)': genbank['PcntIDs'],
        }

        result.append(row)

    columns = [
        'Authors',
        'Title',
        'Journal',
        'Year',
        'PMID',
        'Viruses (AI)',
        'Viruses (GB)',

        'NumSeqs (AI)',
        'NumSeqs (GB)',

        'Hosts (AI)',
        'Hosts (GB)',

        'Specimen (AI)',
        'Specimen (GB)',

        'SampleYr (AI)',
        'SampleYr (GB)',

        'Countries (AI)',
        'Countries (GB)',

        # 'Genes (AI)':
        # 'Genes (GB)':

        'SeqMethod (AI)',
        'SeqMethod (GB)',

        'CloneMethod (AI)',
        'CloneMethod (GB)',

        'GenBank (AI)',
        'GenBank (GB)',

        'AlignLens (GB)',
        'PcntIDs (GB)',
    ]

    for i in result:
        for c in columns:
            if c not in i:
                i[c] = ''

    return pd.DataFrame(result), columns


if __name__ == '__main__':
    main()
