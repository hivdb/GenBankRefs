import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from datetime import datetime
from Utilities import get_logger
import subprocess
import pandas as pd


timestamp = datetime.now().strftime('%m_%d')


class Virus:

    _viruses = {}

    def __init__(self, name):
        self.name = name

        if name in self.__class__._viruses:
            return
        else:
            self.__class__._viruses[name] = self

    @classmethod
    def get_virus(cls, name):
        if name in cls._viruses:
            return cls._viruses[name]
        else:
            return cls._viruses['default']

    @property
    def GENES(self):
        return []

    @property
    def reference_folder(self):
        return Path("ReferenceData") / f"{self.name}"

    @property
    def genbank_file(self):
        return self.reference_folder / f"{self.name}.gb"

    @property
    def BLAST_NA_DB_PATH(self):
        return self.reference_folder / "blast" / f"{self.name}_NA_db"

    @property
    def BLAST_AA_DB_PATH(self):
        return self.reference_folder / "blast" / f"{self.name}_AA_db"

    @property
    def DB_FILE(self):
        return self.output_dir / f"{self.name}.db"

    @property
    def output_dir(self):
        return Path("OutputData") / f"{self.name}"

    @property
    def output_excel_dir(self):
        (self.output_dir / 'excels').mkdir(exist_ok=True)
        return self.output_dir / 'excels'

    @property
    def timestamp_dir(self):
        d = self.output_excel_dir / timestamp
        d.mkdir(exist_ok=True)
        return d

    @property
    def exclude_seq_file(self):
        return self.timestamp_dir / f"{self.name}_Excluded_Seqs_{timestamp}.xlsx"

    @property
    def genbank_raw_ref_file(self):
        return self.timestamp_dir / f"{self.name}_Raw_Ref_{timestamp}.xlsx"

    @property
    def genbank_ref_file(self):
        return self.timestamp_dir / f"{self.name}_Ref_{timestamp}.xlsx"

    @property
    def merged_ref_file(self):
        return self.timestamp_dir / f"{self.name}_Merged_Ref_{timestamp}.xlsx"

    @property
    def genbank_feature_file(self):
        return self.timestamp_dir / f"{self.name}_GenBankFeatures_{timestamp}.xlsx"

    @property
    def genbank_feature_check_file(self):
        return self.timestamp_dir / f"{self.name}__GenBankFeatures_check_{timestamp}.xlsx"

    @property
    def genbank_feature_filled_file(self):
        return self.timestamp_dir / f"{self.name}__GenBankFeatures_filled_{timestamp}.xlsx"

    @property
    def genbank_gene_file(self):
        return self.output_dir / f"{self.name}__GenBankGenes.xlsx"

    @property
    def combined_file(self):
        return self.timestamp_dir / f"{self.name}_Combined_{timestamp}.xlsx"

    @property
    def comparison_file(self):
        return None

    @property
    def pubmed_folder(self):
        return Path("Pubmed") / f"{self.name}"

    @property
    def pubmed_file(self):
        return None

    @property
    def pubmed_additional_from_gb(self):
        return None

    @property
    def pubmed_with_index(self):
        return self.timestamp_dir / f'Pubmed_with_index_{timestamp}.xlsx'

    @property
    def pubmed_genbank_combined(self):
        return self.timestamp_dir / f"{self.name}_P_G_Combined_{timestamp}.xlsx"

    @property
    def genbank_unmatch_file(self):
        return self.timestamp_dir / f"{self.name}_genbank_unmatch_{timestamp}.xlsx"

    @property
    def chord_diagram_file(self):
        return self.timestamp_dir / f"{self.name}_chord1_{timestamp}.html"

    @property
    def chord_diagram_file2(self):
        return self.timestamp_dir / f"{self.name}_chord2_{timestamp}.html"

    @property
    def chord_diagram_file3(self):
        return self.timestamp_dir / f"{self.name}_chord3_{timestamp}.html"

    @property
    def chord_table_file(self):
        return self.timestamp_dir / f"{self.name}_chord_table_{timestamp}.xlsx"

    def get_logger(self, logger_name):
        logger_file = f'{self.name}_datalog_{logger_name}.txt'
        if not getattr(self, logger_file, None):
            setattr(self, logger_file, get_logger(
                self.output_dir / logger_file))
        return getattr(self, logger_file)

    def build_blast_db(self):
        build_blast_db(self)

    def process_features(self, features_df, genes):
        features_df = self._process_features(features_df)

        for i, row in features_df.iterrows():
            g_list = genes[genes['Accession'] == row['Accession']]
            g_list = set(g_list['Gene'].tolist())
            if (g_list == set(self.GENES)):
                features_df.loc[i, 'Genome'] = 1

            g_list = ', '.join(sorted(list(g_list)))
            features_df.loc[i, 'Genes'] = g_list

        features_df['Country'] = features_df[
            'country_region'].str.split(":").str[0]

        features_df = add_feature_from_non_pubmed_paper(self, features_df)

        for i, row in features_df.iterrows():
            if 'patent' in row['Description'].lower():
                features_df.at[i, 'Comment'] = 'patent'
            if 'Modified Microbial Nucleic Acid' in row['Description']:
                features_df.at[i, 'Comment'] = 'Modified Microbial Nucleic Acid'
            if 'CONSTRUCT' in row['Description'].upper():
                features_df.at[i, 'Comment'] = 'CONSTRUCT'
            if 'COMPOSITIONS' in row['Description'].upper():
                features_df.at[i, 'Comment'] = 'CONSTRUCT'
            if 'monoclonal antibody' in row['Description'].lower():
                features_df.at[i, 'Comment'] = 'antibody'
            if 'MICROARRAY' in row['Description'].upper():
                features_df.at[i, 'Comment'] = 'MICROARRAY'
            if 'conformation' in row['Description'].lower():
                features_df.at[i, 'Comment'] = 'conformation'

        return features_df

    def _process_features(self, features_df):
        return features_df

    def process_gene_list(self, gene_df):
        return gene_df


Virus('default')


def build_blast_db(virus):

    reference_aa_file = virus.reference_folder / f"{virus.name}_RefAAs.fasta"

    subprocess.run(
        f"makeblastdb -in {reference_aa_file} "
        f"-dbtype prot -out {virus.BLAST_AA_DB_PATH}",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=True
    )


def add_feature_from_non_pubmed_paper(virus, features_df):
    csv_data = []
    for i in virus.pubmed_folder.iterdir():
        if i.suffix == '.csv':
            csv_data.append(pd.read_csv(i))

    if not csv_data:
        return features_df

    csv_data = pd.concat(csv_data)

    for i, row in features_df.iterrows():
        match = csv_data[csv_data['Accession'] == row['Accession']]
        if not match.empty:
            features_df.at[i, 'Country'] = match['Country'].tolist()[0]
            features_df.at[i, 'collection_date'] = match['IsolateYear'].tolist()[0]
            features_df.at[i, 'isolate_source'] = match['isolate_source'].tolist()[0]

    return features_df
