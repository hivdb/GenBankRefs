from Bio import Entrez
import pandas as pd
import time

Entrez.email = "rshafer.stanford.edu"


def search_pubmed(query):
    time.sleep(0.1)
    handle = Entrez.esearch(db="pubmed", term=query, retmax=3)
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


def search_submission_sets_without_PMID(virus, genbank_unmatched, overwrite=False):

    cache_file = virus.output_excel_dir / f'{virus.name}_pubmed_search.xlsx'

    answers = []

    if cache_file.exists() and not overwrite:
        answers = pd.read_excel(cache_file)
        answer_map = {
            int(i['RefID']): i['PMID']
            for _, i in answers.iterrows()
        }
        for idx, row in genbank_unmatched.iterrows():
            genbank_unmatched.at[idx, 'PMID'] = answer_map.get(row['RefID'], '')

    for idx, row in genbank_unmatched.iterrows():
        if 'PMID' in row and row['PMID']:
            continue

        authors = row['Authors']
        author_pmid = search_pubmed(authors)

        accession_pmids = []
        for i in row['accession'].split(','):
            accession_pmids.extend(search_pubmed(i.strip()))

        pmid = list(sorted(set(author_pmid) | set(accession_pmids)))

        genbank_unmatched.loc[idx, 'PMID'] = ', '.join(
            [str(p) for p in pmid])

        answers.append({
            'RefID': row['RefID'],
            'Title': row['Title'],
            'Authors': row['Authors'],
            'Journal': row['Journal'] if row['Journal'].lower() != 'unpublished' else '',
            'Year': row['Year'] if row['Year'] else '',
            'Accession': row['accession'],
            'PMID': ', '.join([str(p) for p in pmid]),
        })
        print('pubmed search', idx)

    pd.DataFrame(answers).to_excel(cache_file, index=False)
