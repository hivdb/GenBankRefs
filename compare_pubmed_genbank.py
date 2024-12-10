import pandas as pd
from combine.match_pm_gb import match_pm_gb
from combine.combine_file import combine_file, format_table
from combine.sum_genbank import summarize_genbank_by_seq
from combine.sum_genbank import summarize_genbank_by_ref
from combine.sum_pubmed import summarize_pubmed_data
from combine.sum_genbank import summarize_genbank_full_genome
from combine.sum_combined import summarize_combined_data
from combine.translate_value import categorize_host_specimen


def compare_pubmed_genbank(virus_obj):

    if not virus_obj.pubmed_folder.exists():
        print('Pubmed file not found')
        return

    summrize = input('Summarize tables? [y/n]')
    summrize = summrize == 'y'

    genbank_ref_file = virus_obj.combined_file
    genbank_ref = pd.read_excel(genbank_ref_file, dtype=str).fillna('')

    if summrize:
        summarize_genbank_by_ref(genbank_ref)

    genbank_file = virus_obj.genbank_feature_check_file
    genbank = pd.read_excel(genbank_file, dtype=str).fillna('')

    if summrize:
        summarize_genbank_by_seq(genbank)
        summarize_genbank_full_genome(genbank_ref)

    pubmed_file = virus_obj.pubmed_file
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

    if summrize:
        summarize_pubmed_data(pubmed)

    pubmed = pd.concat([pubmed, pd.read_excel(
        virus_obj.pubmed_additional_from_gb, dtype=str).fillna('')],
        ignore_index=True)

    print('#Pubmed Ref with additional GenBank PMID', len(pubmed))
    pubmed['RefID'] = pubmed.index + 1
    pubmed.to_excel(str(virus_obj.pubmed_folder / 'Pubmed.xlsx'), index=False)

    pubmed_match, genbank_match, pubmed_unmatch, genbank_unmatch = match_pm_gb(
        pubmed, genbank_ref)

    # pubmed_unmatch.to_excel('pubmed_unmatch.xlsx')
    # pd.DataFrame(genbank_unmatch).to_excel('genbank_unmatch.xlsx')

    print('Pumbed only', len(pubmed_unmatch))

    if summrize:
        summarize_pubmed_data(pubmed_unmatch)

    print('GenBank only', len(genbank_unmatch))

    if summrize:
        summarize_genbank_by_ref(genbank_unmatch)

    combined, columns = combine_file(
        pubmed_match, pubmed_unmatch, genbank_unmatch)

    combined = categorize_host_specimen(
        combined, 'Hosts (PM)', 'Specimen (PM)')
    combined['Hosts (PM)'] = combined['CleanedHost']
    combined['Specimen (PM)'] = combined['CleanedSpecimen']

    combined.to_excel(str(virus_obj.pubmed_genbank_combined),
                      index=False, columns=columns)
    format_table(str(virus_obj.pubmed_genbank_combined))

    if summrize:
        summarize_combined_data(combined, genbank)
        summarize_genbank_full_genome(genbank_match)

    return pubmed, pubmed_match


if __name__ == '__main__':
    compare_pubmed_genbank()
