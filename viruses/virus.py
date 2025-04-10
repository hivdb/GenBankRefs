from pathlib import Path
from datetime import datetime
from Utilities import get_logger
import subprocess
import pandas as pd
from bioinfo import dump_fasta
from bioinfo import load_fasta
from collections import defaultdict
import matplotlib as mpl
from distinctipy import get_colors
from operator import itemgetter
import matplotlib.pyplot as plt
import numpy as np


timestamp = datetime.now().strftime('%m_%d')


class Virus:

    _viruses = {}

    def __init__(self, name, full_name=None):
        self.name = name
        self.full_name = full_name

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
    def special_accessions(self):
        return []

    @property
    def GENES(self):
        return []

    @property
    def reference_folder(self):
        return Path("ReferenceData") / f"{self.name}"

    @property
    def ref_na_path(self):
        return Path("ReferenceData") / self.name / f"{self.name}_RefNAs.fasta"

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
        # timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        return self.output_dir / f"{self.name}-latest.db"

    @property
    def output_dir(self):
        return Path("OutputData") / f"{self.name}"

    @property
    def alignment_folder(self):
        (self.output_excel_dir / 'alignment').mkdir(exist_ok=True)
        return self.output_excel_dir / 'alignment'

    @property
    def output_excel_dir(self):
        (self.output_dir / 'excels').mkdir(exist_ok=True)
        return self.output_dir / 'excels'

    @property
    def isolate_file(self):
        return self.output_excel_dir / f"{self.name}_isolates.xlsx"

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
    def fixed_ref_id_file(self):
        return self.output_dir / f"{self.name}_fixed_RefID.csv"

    @property
    def fixed_pub_id_file(self):
        return self.output_dir / f"{self.name}_fixed_PubID.csv"

    @property
    def paired_pub_id_ref_id_track(self):
        return self.output_dir / f"{self.name}_paired_PubID_RefID_track.csv"

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
    def genbank_gene_filled_file(self):
        return self.timestamp_dir / f"{self.name}__GenBankGenes_filled_{timestamp}.xlsx"

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
    def pubmed_genbank_hardlink(self):
        return None

    @property
    def pubmed_search_missing(self):
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
    def pubmed_unmatch_file(self):
        return self.timestamp_dir / f"{self.name}_pubmed_unmatch_{timestamp}.xlsx"

    @property
    def pubmed_search_result(self):
        return self.pubmed_folder / f"{self.name}_pubmed_search_checked.xlsx"

    @property
    def AI_search_result(self):
        return self.pubmed_folder / f"{self.name}_AI_search_checked.xlsx"

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
    def chord_table_file1(self):
        return self.timestamp_dir / f"{self.name}_chord_table1_{timestamp}.xlsx"

    @property
    def chord_table_file2(self):
        return self.timestamp_dir / f"{self.name}_chord_table2_{timestamp}.xlsx"

    @property
    def db_dump_folder(self):
        (self.output_dir / 'db_json').mkdir(exist_ok=True)
        return self.output_dir / 'db_json'

    def get_logger(self, logger_name):
        logger_file = f'{self.name}_datalog_{logger_name}.txt'
        if not getattr(self, logger_file, None):
            setattr(self, logger_file,
                    get_logger(self.output_dir / logger_file))
        return getattr(self, logger_file)

    def build_blast_db(self):
        build_blast_db(self)

    def process_features(self, features_df, genes):
        features_df = self._process_features(features_df)

        # add genes for feature
        for i, row in features_df.iterrows():
            g_list = genes[genes['Accession'] == row['Accession']]
            g_list = set(g_list['Gene'].tolist())
            if (g_list == set(self.GENES)):
                features_df.loc[i, 'Genome'] = 1

            g_list = ', '.join(sorted(list(g_list)))
            features_df.loc[i, 'Genes'] = g_list

        country_mapping = {
            "Republic of the Congo": "Democratic Republic of the Congo",
            "United States of America": "United States",
            "USA": "United States",
            "UK": "United Kingdom",
            "Great Britain": "United Kingdom",
            "South Korea": "Korea",
            "North Korea": "Korea",
            "Cote d'Ivoire": "Côte d'Ivoire"
            # "Russia": "Russian Federation",
        }
        features_df['Country'] = features_df['country_region'].str.split(
            ":").str[0]
        features_df['Country'] = features_df['Country'].str.strip()
        features_df['Country'] = features_df['Country'].replace(
            country_mapping)

        features_df = add_feature_from_non_pubmed_paper(self, features_df)

        features_df['isolate_source'] = features_df[
            'isolate_source'].str.capitalize()

        features_df['NonClinical'] = ''
        # for i, row in features_df.iterrows():
        #     if 'patent' in row['Description'].lower():
        #         features_df.at[i, 'NonClinical'] = 'patent'
        #     if 'FDA' in row['Description'].upper():
        #         features_df.at[i, 'NonClinical'] = 'patent'
        #     if 'Modified Microbial Nucleic Acid' in row['Description']:
        #         features_df.at[i, 'NonClinical'] = 'Modified Microbial Nucleic Acid'
        #     if 'CONSTRUCT' in row['Description'].upper():
        #         features_df.at[i, 'NonClinical'] = 'CONSTRUCT'
        #     if 'COMPOSITIONS' in row['Description'].upper():
        #         features_df.at[i, 'NonClinical'] = 'CONSTRUCT'
        #     if 'monoclonal antibody' in row['Description'].lower():
        #         features_df.at[i, 'NonClinical'] = 'antibody'
        #     if 'MICROARRAY' in row['Description'].upper():
        #         features_df.at[i, 'NonClinical'] = 'MICROARRAY'
        #     if 'conformation' in row['Description'].lower():
        #         features_df.at[i, 'NonClinical'] = 'conformation'

        return features_df

    def _process_features(self, features_df):
        return features_df

    def process_gene_list(self, gene_df):
        return gene_df

    def translate_cds_name(self, cds):
        return cds

    @property
    def phylo_folder(self):
        d = self.output_excel_dir / 'phylo'
        d.mkdir(exist_ok=True)
        return d

    @property
    def ref_na_gene_map(self):
        return load_fasta(self.reference_folder / f"{self.name}_RefNAs.fasta")

    @property
    def ref_aa_gene_map(self):
        return load_fasta(self.reference_folder / f"{self.name}_RefNAs.fasta")

    def pick_phylo_sequence(self, genes, picked_genes=[], coverage_pcnt=1):
        return pick_phylo_sequence(self,
                                   genes,
                                   picked_genes,
                                   coverage_pcnt=coverage_pcnt)

    def viz_alignment_coverage(self, gene_df):

        for gene in self.GENES:
            sub_gene_df = gene_df[gene_df['Gene'] == gene]

            ref_na = self.ref_na_gene_map[gene]
            ref_na_length = len(ref_na)

            print(f'Gene {gene}, # Seq', len(sub_gene_df))
            print('# Seq with ins: ', len(sub_gene_df[sub_gene_df['NA_num_ins'] == 1]))
            print('# Seq with del: ', len(sub_gene_df[sub_gene_df['NA_num_del'] == 1]))
            print('# Seq with N: ', len(sub_gene_df[sub_gene_df['num_N'] == 1]))
            print('# Seq with <90% cover: ', len(sub_gene_df[sub_gene_df['NA_length'] < ref_na_length * 0.9]))

            pos_pairs = [
                (row['NA_start'], row['NA_stop'])
                for i, row in sub_gene_df.iterrows()
            ]
            image_folder = self.output_excel_dir / 'alignment_coverage'
            image_folder.mkdir(exist_ok=True)
            image_file_path = image_folder / f'{gene}.png'
            viz_alignment_coverage(image_file_path, gene, pos_pairs)

            image_folder = self.output_excel_dir / 'alignment'
            image_folder.mkdir(exist_ok=True)
            image_file_path = image_folder / f'{gene}_na_length.png'
            viz_histogram(image_file_path, [
                row['NA_length']
                for i, row in sub_gene_df.iterrows()],
                title=f'{gene} NA coverage'
            )

            image_file_path = image_folder / f'{gene}_na_num_N.png'
            viz_histogram(image_file_path, [
                row['num_N']
                for i, row in sub_gene_df.iterrows()],
                title=f"{gene} # N"
            )

            image_file_path = image_folder / f'{gene}_na_num_ins.png'
            viz_histogram(image_file_path, [
                row['NA_num_ins']
                for i, row in sub_gene_df.iterrows()],
                title=f"{gene} # ins"
            )

            image_file_path = image_folder / f'{gene}_na_num_del.png'
            viz_histogram(image_file_path, [
                row['NA_num_del']
                for i, row in sub_gene_df.iterrows()],
                title=f"{gene} # del"
            )


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


def pick_phylo_sequence(virus, genes, picked_genes, coverage_pcnt=1):
    """
    Filters nucleotide sequences based on gene name, clinical relevance, sequence coverage, and uniqueness.
    Extracts metadata (Host, Country, Year) and formats it.
    Saves processed sequences into FASTA and metadata CSV files.
    Optionally runs alignment and phylogenetic tree analysis.
    """
    genes = genes.to_dict(orient='records')

    sampleYr_range = {
        (1900, 1990): '<1990',
        (1991, 2000): '1991-2000',
        (2001, 2010): '2001-2010',
        (2011, 2020): '2011-2020',
        (2021, 2025): '2021-2025',
    }

    num_host = len(set(j['Host'] for j in genes))
    if num_host <= len(mpl.colormaps['Set2'].colors):
        host_palette = mpl.colormaps['Set2'].colors
    else:
        host_palette = get_colors(num_host)

    num_sampleyr = len(set(j['IsolateYear'] for j in genes))
    if num_sampleyr <= len(mpl.colormaps['Set3'].colors):
        sampleyr_palette = mpl.colormaps['Set3'].colors
    else:
        sampleyr_palette = get_colors(num_sampleyr)

    num_country = len(set(j['Country'] for j in genes))
    if num_country <= 20:
        country_palatte = mpl.colormaps['tab20'].colors
    else:
        country_palatte = get_colors(num_country)

    host_color_map = {'NA': '#555555'}
    sampleyr_color_map = {'NA': '#555555'}
    country_color_map = {'NA': '#555555'}

    for gene_name in picked_genes:

        run_command = input(f'Run alignment and phylogenetic tree for {gene_name}? [y/n]')
        if run_command == 'n':
            print('Choose not run phylogenetic tree')
            continue

        ref_na = virus.ref_na_gene_map[gene_name]

        dedup_na = []
        num_dump = 0
        g_list = {}
        metadata = []

        idx = 1
        for j in genes:
            if (j['Gene'] != gene_name):
                continue

            if j['NonClinical']:
                continue

            if int(j['NA_length']):
                if (int(j['NA_length']) < (len(ref_na) * coverage_pcnt)):
                    continue
            elif int(j['NA_raw_length']) < (len(ref_na) * coverage_pcnt * 0.98):
                continue

            if j['NA_raw_seq'] in dedup_na:
                num_dump += 1
                continue

            label = j['Accession']
            # label = f"N{idx}"

            host = j['Host'] if j["Host"] else 'NA'
            host = host.rstrip('*').strip()
            if 'and' in host:
                host = host.split('and', 1)[0].strip()

            country = j['Country'] if j["Country"] else 'NA'
            country = country.rstrip('*').strip()
            if ',' in country:
                country = country.split(',', 1)[0].strip()

            sampleyr = str(j['IsolateYear']) if j["IsolateYear"] else 'NA'
            sampleyr = sampleyr.rstrip('*').strip()
            if '-' in sampleyr:
                sampleyr = sampleyr.split('-')[-1].strip()
            elif '–' in sampleyr:
                sampleyr = sampleyr.split('–')[-1].strip()
            elif ',' in sampleyr:
                sampleyr = sampleyr.split(',')[-1].strip()

            if sampleyr != 'NA':
                sampleyr = int(float(sampleyr))
                for (start, stop), sampleYr_name in sampleYr_range.items():
                    if (sampleyr >= start) and (sampleyr <= stop):
                        sampleyr = sampleYr_name
                        break

            if host == 'NA' or country == 'NA' or sampleyr == 'NA':
                continue

            g_list[label] = j['NA_raw_seq']
            dedup_na.append(j['NA_raw_seq'])

            host_color_map[host] = host_color_map.get(
                host,
                mpl.colors.to_hex(host_palette[len(host_color_map)])
            )
            country_color_map[country] = country_color_map.get(
                country,
                mpl.colors.to_hex(country_palatte[len(country_color_map)])
            )
            sampleyr_color_map[sampleyr] = sampleyr_color_map.get(
                sampleyr,
                mpl.colors.to_hex(sampleyr_palette[len(sampleyr_color_map)])
            )

            metadata.append({
                'label': label,
                'Host': host,
                'Country': country,
                'SampleYr': sampleyr,
                'Host_color': host_color_map[host],
                'Country_color': country_color_map[country],
                'SampleYr_color': sampleyr_color_map[sampleyr],
                # 'Source': j['isolate_source'] if j["isolate_source"] else 'NA',
            })
            idx += 1

        print(f"{virus.name} Gene {gene_name} picked sequence:", len(g_list))
        print(f"{virus.name} Gene {gene_name} unpicked sequence:", len(genes) - len(g_list))
        print(f"{virus.name} Gene {gene_name} duplicated sequence:", num_dump)

        option_one_per_pattern = input('One per pattern? [y/n]') == 'y'

        if option_one_per_pattern:
            metadata, g_list = get_sequences_limited_pattern_each(
                    metadata, g_list, region=True)

        option_country_to_region = input('Change country to region? [y/n]') == 'y'
        if option_country_to_region:
            metadata = convert_country_to_region(metadata)

        print('=' * 80)
        print('# Sequences for building phylogenetic tree:', len(g_list))
        print('=' * 80)

        pd.DataFrame(metadata).to_csv(virus.phylo_folder / f"{gene_name}_metadata.csv", index=False)

        dump_fasta(virus.phylo_folder / f"{gene_name}_ref_na.fasta", {gene_name: ref_na})
        dump_fasta(virus.phylo_folder / f'{gene_name}_isolates.fasta', g_list)

        cmds = (
            f"cd {virus.phylo_folder}; "
            f"rm -rf output_{gene_name}; "
            f"ViralMSA.py -e hivdbteam@list.stanford.edu -s {gene_name}_isolates.fasta -o output_{gene_name} -r {gene_name}_ref_na.fasta --omit_ref; "
            f"cd output_{gene_name}; "
            f"iqtree2 -s {gene_name}_isolates.fasta.aln -m GTR+G4+F -bb 1000 -nt AUTO; "
        )
        print('Run command:')
        print(cmds)

        option_no_output = input("Show tree building details? [y/n]") == 'y'
        if option_no_output:
            subprocess.run(
                cmds,
                # stdout=subprocess.PIPE,
                # stderr=subprocess.PIPE,
                text=True,
                shell=True
            )
        else:
            subprocess.run(
                cmds,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                # text=True,
                shell=True
            )


def get_sequences_limited_pattern_each(metadata, g_list):
    # print(len(g_list))
    pattern_acc = defaultdict(list)

    for i in metadata:

        pattern_acc[(
            i['Host'],
            i['Country'],
            i['SampleYr']
        )].append(i)

    num_pattern = 1
    if len(pattern_acc) < 50:
        num_pattern = 4
    elif (len(pattern_acc) < 100):
        num_pattern = 3
    elif (len(pattern_acc) < 200):
        num_pattern = 2

    metadata = [
        j
        for p, acc_list in pattern_acc.items()
        for j in sorted(acc_list, key=itemgetter('label'))[0: num_pattern]
    ]

    keep_acc = [
        acc['label']
        for acc in metadata
    ]

    g_list = {
        k: v
        for k, v in g_list.items()
        if k in keep_acc
    }

    return metadata, g_list


def convert_country_to_region(metadata):

    # Define subregions and their assigned colors
    subregion_colors = {
        "Central Asia": '#0D4A70',
        "Eastern Asia": '#226E9C',
        "Northern Asia": '#3C93C2',
        "Southern Asia": '#6CB0D6',
        "Western Asia": '#9EC9E2',

        "Eastern Africa": '#06592A',
        "Middle Africa": '#22bb3b',
        "Northern Africa": '#40ad5a',
        "Southern Africa": '#6cba7d',
        "Western Africa": '#9ccea7',

        "Central Europe": '#8f003b',
        "Eastern Europe": '#c40f5b',
        "Northern Europe": '#e32977',
        "Southern Europe": '#e95694',
        "Western Europe": '#ed85b0',


        'Northern America': '#ffffff',
    }
    country_color_map = {'NA': '#555555'}
    from countryinfo import CountryInfo
    mapper = {
        'Yugoslavia': 'Southern Europe',
        'Kosovo': 'Southern Europe',
        'North Macedonia': 'Southern Europe',
    }
    for i in metadata:
        # print(i['Country'])
        if i['Country'] in mapper:
            i['Country'] = mapper[i['Country']]
        else:
            i['Country'] = CountryInfo(i['Country']).subregion()

        country = i['Country']
        country_color_map[country] = subregion_colors[country]
        i['Country_color'] = country_color_map[country]

        i['Country'] = country.split()[-1] + f' ({country.split()[0]})'

    metadata = [
        i
        for i in metadata
        if 'america' not in i['Country'].lower()
    ]

    return metadata


def viz_alignment_coverage(image_file_path, gene, position_pairs):

    position_pairs.sort(key=lambda x: -(x[1] - x[0]))

    plt.figure(figsize=(10, 4))

    for i, (start, end) in enumerate(position_pairs):
        plt.hlines(y=i, xmin=start, xmax=end, color='blue', linewidth=2)

    plt.xlabel('Gene Position (bp)')
    plt.ylabel('# Seq')
    plt.title(f'{gene}')
    plt.tight_layout()

    x = [
        (i[1] - i[0] + 1)
        for i in position_pairs
    ]
    tick_step = closest_smaller_base(max(x) / 5)
    ticks = list(range(0, max(x)+1, tick_step))
    if max(x) not in ticks:
        ticks.append(max(x))
    # print(ticks)
    plt.xticks(ticks, rotation=90)

    plt.savefig(str(image_file_path), dpi=300)
    plt.close()


def closest_smaller_base(n):
    bases = [10, 100, 1000]
    smaller_bases = [b for b in bases if b <= n]
    return max(smaller_bases) if smaller_bases else None


def viz_histogram(image_file_path, data, title):
    data = np.array(data)
    unique_values = np.sort(np.unique(data))
    counts = [np.sum(data == val) for val in unique_values]
    plt.bar(unique_values, counts, edgecolor='black')

    plt.title(title)
    # plt.xlabel('NA length')
    plt.ylabel('# Seq')
    plt.savefig(str(image_file_path), dpi=300)
    plt.close()
