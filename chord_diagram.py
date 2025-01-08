from itertools import combinations
from itertools import product
import pandas as pd


def gen_chord_diagram(virus_obj, combined, features):

    matches = combined[(combined['match'] == 'Yes')]

    accessions = set([
        j.strip()
        for i, v in matches.iterrows()
        for j in v['GenBank (GB)'].split(',')
        ])

    features = features[features['Accession'].isin(list(accessions))]
    # print(len(features))

    columns = [
        'Genes',
        'isolate_source',
        'Host',
        'Country',
        'IsolateYear',
    ]

    year_range = generate_year_ranges(features['IsolateYear'].tolist())
    # print(year_range)
    ordering = get_label_ordering(columns, features, year_range)
    # print(ordering)

    df = []

    host_list = features['Host'].unique()
    gene_list = list(features['Genes'].unique())
    gene_list = [i for i in gene_list if i.strip()]
    gene_list = sorted(list(set([
        j.strip()
        for i in gene_list
        for j in i.split(',')
    ])))
    gene_list.append('')

    for s_value, t_value in product(gene_list, host_list):
        if s_value:
            weight = len(features[
                    (features['Genes'].str.contains(s_value, na=False)) &
                    (features['Host'] == t_value)
                ])
        else:
            weight = len(features[
                (features['Genes'] == s_value) &
                (features['Host'] == t_value)
            ])

        df.append({
                'source': s_value or 'Other Genes',
                'target': t_value or 'Other Host',
                'weight': weight
            })

    columns = [
        'isolate_source',
        'Host',
        'Country',
        'IsolateYear',
    ]

    for s, t in combinations(columns, 2):
        if s == 'isolate_source' and t != 'Host':
            continue

        source_list = features[s].unique()
        target_list = features[t].unique()
        # print(source_list, target_list)
        for s_value, t_value in product(source_list, target_list):
            weight = len(features[
                (features[s] == s_value) &
                (features[t] == t_value)
            ])

            if t == 'IsolateYear' and t_value:
                t_value = get_year_range(t_value, year_range)

            df.append({
                'source': s_value or f'Other {s}',
                'target': t_value or f'Other {t}',
                'weight': weight
            })

    df = pd.DataFrame(df)
    # print(df)

    from d3blocks import D3Blocks
    d3 = D3Blocks()
    d3.chord(df, ordering=ordering, filepath=virus_obj.chord_diagram_file)

    return df


def get_label_ordering(columns, features, year_range):
    ordering = []
    for c in columns:
        column_values = list(features[c].unique())
        if c == 'IsolateYear':
            column_values = list(set([
                get_year_range(i, year_range) if i else i
                for i in column_values
            ]))

        if c == 'Genes':
            other = '' in column_values
            column_values = [
                i for i in column_values if i.strip()
            ]
            column_values = list(set([
                j.strip()
                for i in column_values
                for j in i.split(',')
            ]))
            if other:
                column_values.append('')

        if '' in column_values:
            column_values = sorted([i for i in column_values if i])
            column_values.append(f'Other {c}')

        ordering += column_values

    return ordering


def generate_year_ranges(years, range_length=5):
    if not years:
        return []

    years = [int(y) for y in years if y]

    start_year = min(years)
    stop_year = max(years)

    ranges = []
    current_year = start_year
    while stop_year >= current_year:
        ranges.append((current_year, current_year + range_length - 1))
        current_year += range_length

    return ranges


def get_year_range(year, ranges):
    year = int(year)

    for start, stop in ranges:
        if (start <= year) and (year <= stop):
            return f"{start} - {stop}"

    return None
