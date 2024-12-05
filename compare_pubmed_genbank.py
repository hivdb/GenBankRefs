import pandas as pd
from pathlib import Path
from combine.match_pm_gb import match_pm_gb
from combine.combine_file import combine_file, format_table
from combine.sum_genbank import summarize_genbank_by_seq
from combine.sum_genbank import summarize_genbank_by_ref
from combine.sum_pubmed import summarize_pubmed_data
from combine.sum_genbank import summarize_genbank_full_genome
from combine.sum_combined import summarize_combined_data
from combine.translate_value import categorize_host_specimen


def main():
    folder = Path('./combined')

    genbank_ref_file = folder / 'CCHF_Combined_12_04.xlsx'
    genbank_ref = pd.read_excel(genbank_ref_file, dtype=str).fillna('')
    summarize_genbank_by_ref(genbank_ref)

    genbank_file = folder / 'CCHF__GenBankFeatures_12_04_check.xlsx'
    genbank = pd.read_excel(genbank_file, dtype=str).fillna('')
    summarize_genbank_by_seq(genbank)

    summarize_genbank_full_genome(genbank_ref)

    pubmed_file = folder / 'ReferenceSummary_Dec4.xlsx'
    pubmed = pd.read_excel(pubmed_file, dtype=str).fillna('')
    pubmed = pubmed[
        (pubmed['Resolve Title'] != 'Unlikely') &
        (pubmed['Resolve Seq'] != 'No') &
        (
            (pubmed['Reviewer(s) Seq'] == 'Yes') |
            (pubmed['GPT seq (Y/N)'] == 'Yes')
        )
    ]
    print('#Pubmed Ref', len(pubmed))
    summarize_pubmed_data(pubmed)

    pubmed_additional_ref_file_from_gb = folder / 'ReferenceSummary_Genbank_Nov20.xlsx'
    pubmed = pd.concat([pubmed, pd.read_excel(
        pubmed_additional_ref_file_from_gb, dtype=str).fillna('')],
        ignore_index=True)

    print('#Pubmed Ref', len(pubmed))
    pubmed.to_excel(str(folder / 'Pubmed.xlsx'))

    pubmed_match, genbank_match, pubmed_unmatch, genbank_unmatch = match_pm_gb(
        pubmed, genbank_ref)

    # pubmed_unmatch.to_excel('pubmed_unmatch.xlsx')
    # pd.DataFrame(genbank_unmatch).to_excel('genbank_unmatch.xlsx')

    print('Pumbed only', len(pubmed_unmatch))
    summarize_pubmed_data(pubmed_unmatch)

    print('GenBank only', len(genbank_unmatch))
    summarize_genbank_by_ref(genbank_unmatch)

    combined, columns = combine_file(
        pubmed_match, pubmed_unmatch, genbank_unmatch)

    combined = categorize_host_specimen(combined, 'Hosts (PM)', 'Specimen (PM)')
    combined['Hosts (PM)'] = combined['CleanedHost']
    combined['Specimen (PM)'] = combined['CleanedSpecimen']

    combined.to_excel(str(folder / 'combined.xlsx'), index=False, columns=columns)
    format_table(str(folder / 'combined.xlsx'))

    summarize_combined_data(combined, genbank)

    summarize_genbank_full_genome(genbank_match)


if __name__ == '__main__':
    main()
