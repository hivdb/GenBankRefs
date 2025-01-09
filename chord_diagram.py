from itertools import combinations
from itertools import product
import pandas as pd
from copy import deepcopy


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

    column_color = {
        'Genes': '#2b7d9a',
        'isolate_source': '#4ba6c5',
        'Host': '#c574b0',
        'Country': '#a3ad41',
        'IsolateYear': '#778322',
    }

    year_range = generate_year_ranges(features['IsolateYear'].tolist())
    # print(year_range)
    ordering = get_label_ordering(columns, features, year_range)
    # print(ordering)

    df1 = get_gene_host_link(features, column_color)
    df2 = get_paired_link(features, column_color, ['isolate_source', 'Host'], year_range)

    columns = [
        'Host',
        'Country',
        'IsolateYear',
    ]
    df3 = get_paired_link(features, column_color, columns, year_range)

    total = len(features)
    draw_figure(virus_obj.chord_diagram_file,
                deepcopy(df1) + deepcopy(df2) + deepcopy(df3), total, ordering)
    draw_figure(virus_obj.chord_diagram_file2, df1 + df2, total, ordering)
    draw_figure(virus_obj.chord_diagram_file3, df3, total, ordering)

    get_chord_table(virus_obj.chord_table_file, features, year_range)

    return df1 + df2 + df3


def get_chord_table(save_path, features, year_range):
    columns = ['Host', 'Country', 'IsolateYear']
    df = []
    for (a, b, c) in [columns]:

        list1 = features[a].unique()
        list2 = features[b].unique()
        list3 = features[c].unique()
        # print(source_list, target_list)

        for a1, b1, c1 in product(list1, list2, list3):
            weight = len(features[
                (features[a] == a1) &
                (features[b] == b1) &
                (features[c] == c1)
            ])

            if c1:
                c1 = get_year_range(c1, year_range)

            if not weight:
                continue

            df.append({
                'Host': a1,
                'Country': b1,
                'IsolateYear': c1,
                '#': weight,
                'T': len(features),
                '%': round(weight / len(features) * 100),
            })

    pd.DataFrame(df).to_excel(save_path)


def get_gene_host_link(features, column_color):
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

    for s, t in [('Genes', 'Host')]:
        for s_value, t_value in product(gene_list, host_list):
            if s_value:
                weight = len(features[
                        (features[s].str.contains(s_value, na=False)) &
                        (features[t] == t_value)
                    ])
            else:
                weight = len(features[
                    (features[s] == s_value) &
                    (features[t] == t_value)
                ])

            if not weight:
                continue

            df.append({
                    'source': s_value or f"Other {s}",
                    'target': t_value or f"Other {t}",
                    'weight': weight,
                    'source_column': s,
                    'target_column': t,
                    'source_color': column_color[s],
                    'target_color': column_color[t]
                })

    return df


def get_paired_link(features, column_color, columns, year_range):
    df = []
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

            if not weight:
                continue

            df.append({
                'source': s_value or f"Other {s}",
                'target': t_value or f"Other {t}",
                'weight': weight,
                'source_column': s,
                'target_column': t,
                'source_color': column_color[s],
                'target_color': column_color[t]
            })
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
            other = ('' in column_values)
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

        column_values = sorted([i for i in column_values if i])
        column_values.append(f'Other {c}')

        ordering += column_values

    return ordering


def generate_year_ranges(years, range_length=5, split=3):
    if not years:
        return []

    years = [int(y) for y in years if y]

    start_year = min(years)
    stop_year = max(years)

    if split:
        step = (stop_year - start_year + 1) // split
        ranges = []
        for i in range(split):
            range_start = start_year + i * step

            range_end = range_start + step - 1
            if i == split - 1:
                range_end = stop_year
            ranges.append((range_start, range_end))
    else:

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


def update_other_labels(df, total, pcnt=0):
    other_source = get_other_labels(df, 'source', total, pcnt)
    other_target = get_other_labels(df, 'target', total, pcnt)
    other_list = set(other_source + other_target)
    for i in df:
        if i['source'] in other_list:
            i['source'] = f'Other {i["source_column"]}'
        if i['target'] in other_list:
            i['target'] = f'Other {i["target_column"]}'

        del i['source_column']
        del i['target_column']

    return df


def get_node_color(df):
    node_color = {}
    for i in df:
        node_color[i['source']] = i['source_color']
        node_color[i['target']] = i['target_color']

        del i['source_color']
        del i['target_color']

    return node_color


def get_other_labels(df, label_type, total, pcnt=0):
    other_label_list = []
    labels = set([
        i[label_type]
        for i in df
        if (i[f'{label_type}_column'] != 'IsolateYear')
    ])
    for i in labels:
        if i.startswith('Other'):
            other_label_list.append(i)
            continue

        rows = [
            row['weight']
            for row in df
            if row[label_type] == i
        ]
        if (sum(rows) / total) <= pcnt:
            other_label_list.append(i)

    return other_label_list


def update_ordering(df, ordering):
    sources = [
        row['source']
        for row in df
    ]
    targets = [
        row['target']
        for row in df
    ]

    ordering = [j for j in ordering if j in sources or j in targets]

    return ordering


def draw_figure(save_path, df, total, ordering):
    df = update_other_labels(df, total=total, pcnt=0.1)
    ordering = update_ordering(df, ordering)
    node_color = get_node_color(df)

    df = pd.DataFrame(df)
    # print(df)

    from d3blocks import D3Blocks
    d3 = D3Blocks()
    d3.chord(
        df,
        ordering=ordering,
        arrowhead=-1,
        showfig=False)

    for idx, link in df.iterrows():
        label = d3.node_properties[d3.node_properties['label'] == link['source']]
        d3.node_properties.loc[label.index, 'color'] = node_color[link['source']]

        label = d3.node_properties[d3.node_properties['label'] == link['target']]
        d3.node_properties.loc[label.index, 'color'] = node_color[link['target']]

    d3.set_edge_properties(df, color='source', opacity='source')
    d3.show(filepath=save_path)
