import pandas as pd
from .utils import merge_genbank_list_columns
from openpyxl import load_workbook
from openpyxl.styles import Alignment


def combine_file(pubmed_match, pubmed_unmatch, genbank_unmatch):

    result = []
    for pubmed, genbank_list in pubmed_match:
        row = {
            'Authors': pubmed['Authors'],
            'Title': pubmed['Title'],
            'Journal': pubmed['Journal'],
            'Year': pubmed['Year'],
            'PMID': pubmed['PMID'],

            'Reviewer(s) Seq': pubmed['Reviewer(s) Seq'],
            'GPT seq (Y/N)': pubmed['GPT seq (Y/N)'],
            'Resolve': pubmed['Resolve'],
            'match': 'Yes',

            'Viruses (PM)': pubmed['Viruses'],
            'NumSeqs (PM)': pubmed['NumSeqs'],
            'Hosts (PM)': pubmed['Host'],
            'Specimen (PM)': pubmed['IsolateType'],

            'SampleYr (PM)': pubmed['SampleYr'],
            'Countries (PM)': pubmed['Country'],
            'Genes (PM)': pubmed['Gene'],
            'SeqMethod (PM)': pubmed['SeqMethod'],
            'CloneMethod (PM)': pubmed['CloneMethod'],
            'GenBank (PM)': pubmed['GenBank'],

            'Viruses (GB)': merge_genbank_list_columns(genbank_list, 'Organisms'),
            'NumSeqs (GB)': len([
                i
                for genbank in genbank_list
                for i in genbank['accession'].split(',')
            ]),
            'Hosts (GB)': merge_genbank_list_columns(genbank_list, 'Hosts'),

            'Specimen (GB)': merge_genbank_list_columns(genbank_list, 'Specimens'),
            'SampleYr (GB)':  merge_genbank_list_columns(genbank_list, 'IsolateYears'),
            'Countries (GB)': merge_genbank_list_columns(genbank_list, 'Countries'),
            'Genes (GB)': merge_genbank_list_columns(genbank_list, 'Gene'),
            'SeqMethod (GB)': '',
            'CloneMethod (GB)': '',
            'GenBank (GB)': ', '.join([i.strip().split('.')[0]
                                       for genbank in genbank_list
                                       for i in genbank['accession'].split(',')
                                       ]),
            'AlignLens (GB)': merge_genbank_list_columns(genbank_list, 'AlignLens'),
            'PcntIDs (GB)': merge_genbank_list_columns(genbank_list, 'PcntIDs'),
        }

        result.append(row)

    for r, pubmed in pubmed_unmatch.iterrows():
        row = {
            'Authors': pubmed['Authors'],
            'Title': pubmed['Title'],
            'Journal': pubmed['Journal'],
            'Year': pubmed['Year'],
            'PMID': pubmed['PMID'],

            'Reviewer(s) Seq': pubmed['Reviewer(s) Seq'],
            'GPT seq (Y/N)': pubmed['GPT seq (Y/N)'],
            'Resolve': pubmed['Resolve'],

            'Viruses (PM)': pubmed['Viruses'],
            'NumSeqs (PM)': pubmed['NumSeqs'],
            'Hosts (PM)': pubmed['Host'],
            'Specimen (PM)': pubmed['IsolateType'],
            'SampleYr (PM)': pubmed['SampleYr'],
            'Countries (PM)': pubmed['Country'],
            'Genes (PM)': pubmed['Gene'],
            'SeqMethod (PM)': pubmed['SeqMethod'],
            'CloneMethod (PM)': pubmed['CloneMethod'],
            'GenBank (PM)': pubmed['GenBank'],
        }

        result.append(row)

    for row, genbank in genbank_unmatch.iterrows():
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
            'Specimen (GB)': genbank['Specimens'],
            'SampleYr (GB)': genbank['IsolateYears'],
            'Countries (GB)': genbank['Countries'],
            'Genes (GB)': genbank['Gene'],
            'SeqMethod (GB)': '',
            'CloneMethod (GB)': '',
            'GenBank (GB)': ', '.join([
                i.strip().split('.')[0]
                for i in genbank['accession'].split(',')
            ]),
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

        'Reviewer(s) Seq',
        'GPT seq (Y/N)',
        'Resolve',
        'match',

        'Viruses (PM)',
        'Viruses (GB)',

        'NumSeqs (PM)',
        'NumSeqs (GB)',

        'Hosts (PM)',
        'Hosts (GB)',

        'Specimen (PM)',
        'Specimen (GB)',

        'SampleYr (PM)',
        'SampleYr (GB)',

        'Countries (PM)',
        'Countries (GB)',

        'Genes (PM)',
        'Genes (GB)',

        'SeqMethod (PM)',
        'SeqMethod (GB)',

        'CloneMethod (PM)',
        'CloneMethod (GB)',

        'GenBank (PM)',
        'GenBank (GB)',

        'AlignLens (GB)',
        'PcntIDs (GB)',
    ]

    for i in result:
        for c in columns:
            if c not in i:
                i[c] = ''

    return pd.DataFrame(result), columns


def format_table(excel_file):
    wb = load_workbook(excel_file)
    ws = wb.active

    for row in ws.iter_rows(
            min_row=1, max_row=ws.max_row, min_col=1, max_col=ws.max_column):
        for cell in row:
            cell.alignment = Alignment(
                horizontal="left", vertical="top", wrap_text=True)

    for col in ws.columns:
        col_letter = col[0].column_letter
        ws.column_dimensions[col_letter].width = 20

    wb.save(excel_file)
