import pandas as pd
from combine.match_pm_gb import match_pm_gb
from combine.combine_file import combine_file, format_table
from combine.sum_genbank import summarize_genbank_by_seq
from combine.sum_genbank import summarize_genbank_by_ref
from combine.sum_pubmed import summarize_pubmed
from combine.sum_genbank import summarize_genbank_full_genome
from combine.sum_combined import summarize_combined_data


def compare_pubmed_genbank(virus_obj):

    if not virus_obj.pubmed_folder.exists():
        print('Pubmed file not found')
        return pd.DataFrame(), []

    summrize = input('Summarize tables? [y/n]')
    summrize = summrize == 'y'
    print('\n')

    genbank_ref = pd.read_excel(virus_obj.merged_ref_file, dtype=str).fillna('')

    genbank_feature = pd.read_excel(
        virus_obj.genbank_feature_check_file, dtype=str).fillna('')
    genbank_feature['Genes'] = genbank_feature['Genes'].apply(
        virus_obj.translate_gene)

    genbank_genes = pd.read_excel(
        virus_obj.genbank_gene_file, dtype=str).fillna('')

    if summrize:
        summarize_genbank_by_ref(genbank_ref, virus_obj.genbank_logger)
        summarize_genbank_by_seq(
            genbank_feature, genbank_genes, virus_obj.genbank_logger)
        summarize_genbank_full_genome(
            genbank_ref, genbank_feature, virus_obj.genbank_logger,
            full_gene_set=virus_obj.GENES)

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

    print('Pubmed Literatures:', len(pubmed))
    pubmed = virus_obj.process_pubmed(pubmed)

    if summrize:
        summarize_pubmed(pubmed, virus_obj.pubmed_logger)

    if virus_obj.pubmed_additional_from_gb:
        additional_pubmed = pd.read_excel(
            virus_obj.pubmed_additional_from_gb, dtype=str).fillna('')

        additional_pubmed = virus_obj.process_pubmed(additional_pubmed)
        pubmed = pd.concat([pubmed, additional_pubmed], ignore_index=True)
        print(
            '#Pubmed Literature with additional Literature from GenBank Only', len(pubmed))

    # IF inlude additional PMID from GenBank
    # if summrize:
    #     summarize_pubmed(pubmed, virus_obj.pubmed_logger)

    pubmed['LitID'] = pubmed.index + 1
    pubmed.to_excel(str(virus_obj.pubmed_folder / 'Pubmed.xlsx'), index=False)

    pubmed_match, genbank_match, pubmed_unmatch, genbank_unmatch = match_pm_gb(
        pubmed, genbank_ref, virus_obj.pm_gb_logger)

    # pubmed_unmatch.to_excel('pubmed_unmatch.xlsx')
    # pd.DataFrame(genbank_unmatch).to_excel('genbank_unmatch.xlsx')

    print('Pumbed only Literatures:', len(pubmed_unmatch))

    if summrize:
        virus_obj.pm_gb_logger.info('Pubmed Only:')
        summarize_pubmed(pubmed_unmatch, virus_obj.pm_gb_logger)

    print('GenBank only References:', len(genbank_unmatch))

    if summrize:
        virus_obj.pm_gb_logger.info('GenBank Only:')
        summarize_genbank_by_ref(genbank_unmatch, virus_obj.pm_gb_logger)

    combined, columns = combine_file(
        pubmed_match, pubmed_unmatch, genbank_unmatch,
        genbank_feature, genbank_genes)

    combined.to_excel(str(virus_obj.pubmed_genbank_combined),
                      index=False, columns=columns)
    format_table(str(virus_obj.pubmed_genbank_combined))

    if summrize:
        summarize_combined_data(
            combined, genbank_feature, genbank_genes, virus_obj.pm_gb_logger)
        summarize_genbank_full_genome(
            genbank_match, genbank_feature, virus_obj.pm_gb_logger,
            virus_obj.GENES)

    return pubmed, pubmed_match


if __name__ == '__main__':
    compare_pubmed_genbank()
