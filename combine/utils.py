from collections import defaultdict
import pandas as pd
from collections import Counter


def count_number(rows, key=None):
    if key:
        column_values = [row[key] for row in rows]
    else:
        column_values = rows
    column_values = [
        k if k and pd.notna(k) and pd.notnull(k) else 'NA'
        for k in column_values
    ]
    counter = dict(Counter(column_values))

    return ', '.join([
        f'{k} ({v})'
        for k, v in sorted(counter.items(), key=lambda x: int(x[-1]), reverse=True)
    ])


def merge_genbank_list_columns(genbank_list, key, translater=lambda x: x):

    column_values = [genbank[key] for genbank in genbank_list]
    value_count_list = [
        j.strip()
        for i in column_values
        for j in i.split(',')
    ]
    value_count_list = [
        split_value_count(j)
        for j in value_count_list
    ]

    result = defaultdict(int)
    for value, count in value_count_list:
        result[
            translater(value)
            if value != 'NA' else value
        ] += int(count)

    return ','.join([
        f'{k} ({v})'
        for k, v in sorted(result.items(), key=lambda x: int(x[-1]), reverse=True)
    ])


def split_value_count(value_count):

    value_count = value_count[::-1]

    count = value_count[value_count.find(
        ')') + 1: value_count.find('(')][::-1].strip()
    value = value_count[value_count.find('(') + 1:][::-1].strip()

    return value, count


def split_value_by_comma(df, key):
    country_list = []

    for i, row in df.iterrows():
        country = [i.strip() for i in row[key].split(',')]
        country_list.extend(country)

    return country_list
