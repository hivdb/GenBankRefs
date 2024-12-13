from collections import defaultdict
import pandas as pd
from collections import Counter


def count_rev_sorter(value_count_list):
    return sorted(value_count_list, key=lambda x: int(x[-1]), reverse=True)


def alphabetical_sorter(value_count_list):
    return sorted(value_count_list, key=lambda x: x[0])


def int_sorter(value_count_list):
    return sorted(value_count_list, key=lambda x: int(x[0]) if (x[0] and x[0] != 'NA') else -1)


def count_number(rows, key=None, translater=lambda x: x, sorter=count_rev_sorter):
    if key:
        column_values = [row[key] for row in rows]
    else:
        column_values = rows
    column_values = [
        k if k and pd.notna(k) and pd.notnull(k) else 'NA'
        for k in column_values
    ]
    column_values = [translater(x) for x in column_values]
    counter = dict(Counter(column_values))

    return '\n'.join([
        f'{k} ({v})'
        for k, v in sorter(counter.items())
    ])


def split_value_count(value_count):

    value_count = value_count[::-1]

    count = value_count[value_count.find(
        ')') + 1: value_count.find('(')][::-1].strip()
    value = value_count[value_count.find('(') + 1:][::-1].strip()

    return value, count


def get_values_of_value_count_list(value_count_str):
    count_list = []
    value_list = []
    for i in value_count_str.split(','):
        value, count = split_value_count(i)
        count_list.append(int(count))
        value_list.extend([value] * int(count))

    return value_list


def split_value_by_comma(df, key):
    country_list = []

    for i, row in df.iterrows():
        country = [i.strip() for i in row[key].split(',')]
        country_list.extend(country)

    return country_list


def sum_value_count(value_count_str):
    count_list = []
    value_list = []
    for i in value_count_str.split(','):
        value, count = split_value_count(i)
        count_list.append(int(count))
        value_list.extend([value] * int(count))

    print(sum(count_list))
    # print(sum([int(i) for i in value_list]))
