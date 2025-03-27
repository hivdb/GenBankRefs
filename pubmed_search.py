from Bio import Entrez
import pandas as pd
import time
from itertools import combinations

Entrez.email = "rshafer.stanford.edu"


def search_pubmed(query, db='pubmed', retmax=1):
    time.sleep(1)
    handle = Entrez.esearch(db=db, term=query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]


def fetch_pubmed_details(pubmed_ids):
    handle = Entrez.efetch(
        db="pubmed", id=",".join(pubmed_ids),
        rettype="abstract",
        retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    articles = []
    for article in records["PubmedArticle"]:
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        # abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No abstract"])[0]
        authors = [
            author["LastName"] + " " + author["ForeName"]
            for author in article["MedlineCitation"]["Article"].get("AuthorList", [])
            if "LastName" in author and "ForeName" in author]
        pubmed_id = article["MedlineCitation"]["PMID"]
        articles.append({
            "PubMed ID": pubmed_id,
            "Title": title,
            "Authors": ", ".join(authors),
            # "Abstract": abstract
        })

    return articles


def search_by_pubmed_API(
        virus, genbank, overwrite=False):

    if virus.pubmed_search_result.exists() and not overwrite:
        return pd.read_excel(virus.pubmed_search_result)

    cache_file = virus.output_excel_dir / f'{virus.name}_pubmed_search.xlsx'

    answers = []

    if cache_file.exists():
        answers = pd.read_excel(cache_file)
        answer_map = {
            int(i['RefID']): i['PMID']
            for _, i in answers.iterrows()
        }
        for idx, row in genbank.iterrows():
            genbank.at[idx, 'PMID'] = answer_map.get(row['RefID'], '')

        answers = answers.to_dict(orient='records')

    genbank['PMID'] = genbank['PMID'].fillna('').astype(str)

    for idx, row in genbank.iterrows():
        if 'PMID' in row and row['PMID'].strip():
            continue

        title_pmid = []
        title = row['Title'].replace('Direct Submission', '').strip()
        if title:
            title_pmid = search_pubmed(title, retmax=1)

        author_pmid = []
        authors = row['Authors']
        if authors != 'NCBI':
            author_pmid += search_pubmed(authors, retmax=5)
            authors = list(set([
                i.strip()
                for i in authors.split(',')
            ]))
            combo = list(combinations(authors, 2))
            # print('# Combo', len(combo))
            for a1, a2 in combo:
                search_term = f'{a1} and {a2} and {virus.name}'
                author_pmid += search_pubmed(search_term, retmax=5)

        accession_pmids = []
        for i in row['accession'].split(','):
            accession_pmids.extend(search_pubmed(i.strip(), retmax=3))

        # accession_pmids_2 = []
        # for i in row['accession'].split(','):
        #     accession_pmids_2.extend(search_pubmed(i.strip(), db='pmc'))

        pmid = list(sorted(set(author_pmid) | set(accession_pmids)))

        genbank.loc[idx, 'PMID'] = ', '.join(
            [str(p) for p in pmid])

        answers.append({
            'RefID': row['RefID'],
            'Title': row['Title'],
            'Authors': row['Authors'],
            'Journal': row['Journal'] if row['Journal'].lower() != 'unpublished' else '',
            'Year': row['Year'] if row['Year'] else '',
            'Accession': row['accession'],
            'PMID': ', '.join([str(p) for p in pmid]),
            'PMID_title': ', '.join([str(p) for p in title_pmid]),
            'PMID_author': ', '.join([str(p) for p in set(author_pmid)]),
            'PMID_acc': ', '.join([str(p) for p in set(accession_pmids)]),
            # 'PMID_acc2': ', '.join([str(p) for p in set(accession_pmids_2)]),
        })
        print('pubmed search', row['RefID'])

        pd.DataFrame(answers).to_excel(cache_file, index=False)

    print('Please check PubMed Search result by hand.')
    exit()
