import pandas as pd
from pathlib import Path
from combine.match_pm_gb import match_pm_gb
from combine.combine_file import combine_file, format_table
from combine.sum_genbank_by_seq import summarize_genbank_by_seq
from combine.sum_pubmed import summarize_pubmed_data
from combine.sum_combined import summarize_combined_data

def main():
    folder = Path('./combined')
    genbank_file = folder / 'CCHF_Combined_11_21.xlsx'
    pubmed_file = folder / 'ReferenceSummary_Nov25.xlsx'
    genbank_ref_file = folder / 'ReferenceSummary_Genbank_Nov20.xlsx'

    genbank = pd.read_excel(genbank_file, dtype=str).fillna('')
    pubmed = pd.read_excel(pubmed_file, dtype=str).fillna('')

    pubmed = pubmed[((pubmed['Reviewer(s) Seq'] == 'Yes') & (
        pubmed['GPT seq (Y/N)'] == 'Yes')) | (pubmed['Resolve'] == 'Yes')]

    # pubmed = pd.concat([pubmed, pd.read_excel(genbank_ref_file, dtype=str).fillna('')])
    # print(len(pubmed))
    # raise

    pubmed.to_excel(str(folder / 'Pubmed.xlsx'))

    pubmed_match, pubmed_unmatch, genbank_unmatch = match_pm_gb(pubmed, genbank)

    # pubmed_unmatch.to_excel('pubmed_unmatch.xlsx')
    # pd.DataFrame(genbank_unmatch).to_excel('genbank_unmatch.xlsx')

    combined, columns = combine_file(pubmed_match, pubmed_unmatch, genbank_unmatch)
    combined.to_excel(str(folder / 'combined.xlsx'), index=False, columns=columns)
    format_table(str(folder / 'combined.xlsx'))

    print('Combined')
    summarize_combined_data(combined)

    print("PubMed")
    summarize_pubmed_data(pubmed)

    print("GenBank")
    summarize_genbank_by_seq(genbank)


if __name__ == '__main__':
    main()
