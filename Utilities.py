import pandas as pd
import re
from collections import Counter
import logging
from statistics import median
from collections import defaultdict
from pathlib import Path
import csv


def get_logger(logging_file):
    logger = logging.getLogger('logger')
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(message)s')

    file_handler = logging.FileHandler(logging_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # console_handler = logging.StreamHandler()
    # console_handler.setLevel(logging.DEBUG)
    # console_handler.setFormatter(formatter)
    # logger.addHandler(console_handler)

    fd = open(logging_file, 'w')

    class SimpleLogger:

        def info(self, *args, **kargs):
            print(*args, **kargs, file=fd)

        def report(self, report):
            for section in report:
                for pid, part in enumerate(section):
                    if isinstance(part, tuple):
                        self.info(*part)
                    elif isinstance(part, list):
                        self.info(*part)
                    else:
                        self.info(part)
                    if pid < len(section) - 1:
                        self.info('-' * 80)
                self.info('=' * 80)

    return SimpleLogger()


def process_author_field(names):
    if names == '':
        return 'NCBI'

    names = names.replace(" and", ",")
    name_list = names.split(".,")

    processed_names = []
    for name in name_list:
        parts = name.strip().split(",", 1)

        last_name = parts[0]  # The last name is the first part
        last_name = last_name.capitalize()
        initials = parts[1]  # The initials are the second part
        # Take only the first character of the initials
        first_initial = initials[0]
        processed_name = last_name + ' ' + first_initial + '.'
        # processed_name = last_name + ' ' + initials.replace('.', '')
        processed_names.append(processed_name)

    processed_names = ', '.join(processed_names)
    return (processed_names)


def extract_year_from_journal(text):
    # Match the last (xxxx), because in the middle it can mean page number
    matches = re.findall(r'\((\d{4})\)', text)
    if matches:
        return matches[-1]

    match = re.search(r'\d{2}-[A-Z]{3}-(\d{4})', text)
    if match:
        return match.group(1)
    return ''


def extract_year_from_date_fields(text):
    text = str(text) if text is not None else ""
    match = re.search(r"\d{4}", text)
    if match:
        return match.group()
    else:
        return ""


# Each list is a string in which items are separated by ', '
def combine_items_in_different_lists(lists, spliter=None):
    unique_value = set()

    for value in lists:
        value = str(value)
        if spliter:
            value_list = [
                i.strip()
                for i in value.split(spliter)
                if i.strip()
            ]
        else:
            value_list = [value]

        value_list = [
            v
            for v in value_list
            if v and not any([v in u for u in unique_value])
        ]
        for v in value_list:
            unique_value.add(v)

    unique_items = ', '.join(sorted(list(unique_value)))
    return unique_items


def combine_authors_in_different_lists(lists, spliter=None):

    lists = sorted(lists)

    firsts = []
    middles = []
    lasts = []
    for value in lists:
        value = str(value)
        value_list = [
                i.strip()
                for i in value.split(spliter)
                if i.strip()
            ]
        if len(value_list) < 1:
            first = []
            middle = []
            last = []
        elif len(value_list) == 1:
            first = value_list
            middle = []
            last = []
        else:
            first = value_list[:1]
            middle = value_list[1:-1]
            last = value_list[-1:]

        firsts.extend(first)
        middles.extend(middle)
        lasts.extend(last)

    firsts = list(dict.fromkeys(firsts))
    lasts = list(dict.fromkeys(lasts))
    middles = [
        i
        for i in middles
        if (i not in firsts) and (i not in lasts)
    ]
    # middles = list(dict.fromkeys(middles))

    authors = firsts + middles + lasts
    return ', '.join(authors)


# Convert {key: [values], key: [values]] to a list of sets comprising {key, values}
# In this program, each key is the index of a row and the values contain
# one or more indexes of rows that share some property
# The subset flag code ensures sets do not share items with one another

# TODO: Algorithm
def convert_dict_to_list_of_sets(df, matched_indexes):
    list_of_sets = [
        set([i] + j_list)
        for i, j_list in matched_indexes.items()
    ]
    # print(list_of_sets)

    def get_linked_pair(components):
        for idx, iline in enumerate(components):
            for jdx, jline in enumerate(components):
                if idx >= jdx:
                    continue
                if iline & jline:
                    return idx, jdx

    def split_row_by_pmid(row):
        row_with_pmid = defaultdict(list)
        for i, r in row.iterrows():
            pmid = r['PMID'].strip()
            if not pmid:
                continue
            row_with_pmid[pmid].append(r['RowID'])

        row_wo_pmid = [
            r['RowID']
            for i, r in row.iterrows()
            if not r['PMID'].strip()
        ]

        new_arrays = []
        for i, (pmid, r_list) in enumerate(row_with_pmid.items()):
            if i == 0:
                new_arrays.append(set(r_list + row_wo_pmid))
            else:
                new_arrays.append(set(r_list))

        if not new_arrays:
            new_arrays = [set(row_wo_pmid)]

        return new_arrays

    nlist_of_sets = []
    for iline in list_of_sets:
        row_i = df[df['RowID'].isin(list(iline))]
        row_i_arrays = split_row_by_pmid(row_i)
        nlist_of_sets.extend(row_i_arrays)

    list_of_sets = nlist_of_sets

    pair = get_linked_pair(list_of_sets)
    while pair:
        idx, jdx = pair

        iline = list_of_sets[idx]
        jline = list_of_sets[jdx]

        row_i = df[df['RowID'].isin(list(iline))]
        row_i_arrays = split_row_by_pmid(row_i)

        if len(row_i_arrays) > 1:
            list_of_sets[idx] = row_i_arrays[0]
            list_of_sets = list_of_sets[:idx] + row_i_arrays[1:] + list_of_sets[idx:]
            continue

        row_j = df[df['RowID'].isin(list(jline))]
        row_j_arrays = split_row_by_pmid(row_j)
        if len(row_j_arrays) > 1:
            list_of_sets[jdx] = row_j_arrays[0]
            list_of_sets = list_of_sets[:jdx] + row_j_arrays[1:] + list_of_sets[jdx:]
            continue

        iline = row_i_arrays[0]
        jline = row_j_arrays[0]

        PMID_i = [pi for pi in df[df['RowID'].isin(list(iline))]['PMID'].to_list() if pi]
        PMID_j = [pj for pj in df[df['RowID'].isin(list(jline))]['PMID'].to_list() if pj]

        if not (PMID_i and PMID_j):
            list_of_sets = merge_two_row(list_of_sets, idx, jdx)
        else:
            PMID_i = PMID_i[0]
            PMID_j = PMID_j[0]
            if PMID_i == PMID_j:
                list_of_sets = merge_two_row(list_of_sets, idx, jdx)
            else:
                list_of_sets[jdx] = jline - iline

        pair = get_linked_pair(list_of_sets)

    list_of_rows = [
        j
        for i in list_of_sets
        for j in i
    ]

    return (list_of_sets, list_of_rows)


def merge_two_row(alist, i, j):

    a = alist[i]
    b = alist[j]

    # Remove elements at original indices without shifting
    new_list = [v for idx, v in enumerate(alist) if idx not in (i, j)]
    new_list.append(a | b)
    return new_list


def get_pcnt_authors_overlap(authors1, authors2):
    if len(authors1) == 0 or len(authors2) == 0:
        return 0
    set1 = set(authors1.split(', '))
    set2 = set(authors2.split(', '))

    if (set1 in set2) or (set2 in set1):
        return 1
    shared_set = set1 & set2

    # combined_set = set1 | set2
    # pcnt = len(shared_set) / len(combined_set)

    return max(len(shared_set) / len(set1), len(shared_set) / len(set2))


def get_pcnt_shared_accessions(list1, list2):
    if len(list1) == 0 or len(list2) == 0:
        return 0
    pcnt_shared_accessions = len(set(list1) & set(
        list2)) / len(set(list1) | set(list2))
    return pcnt_shared_accessions


def get_pcnt_shared_stems(list1, list2, stem_length):
    if len(list1) == 0 or len(list2) == 0:
        return 0
    acc_num_stem_list1 = []
    acc_num_stem_list2 = []
    for acc_num in list1:
        acc_num_stem_list1.append(acc_num[:stem_length])
    for acc_num in list2:
        acc_num_stem_list2.append(acc_num[:stem_length])
    set_acc_num_stem1 = set(acc_num_stem_list1)
    set_acc_num_stem2 = set(acc_num_stem_list2)
    pcnt_shared_stems = (
        len(set_acc_num_stem1 & set_acc_num_stem2) /
        len((set_acc_num_stem1 | set_acc_num_stem2))
    )
    return pcnt_shared_stems


def calc_year_dif(list1, list2):
    years1 = list1.split(", ")
    years2 = list2.split(", ")
    years1 = [int(year) for year in years1 if year]
    years2 = [int(year) for year in years2 if year]
    max_dif = 0
    for year1 in years1:
        for year2 in years2:
            if abs(year1 - year2) > max_dif:
                max_dif = abs(year1 - year2)
    return max_dif


def count_unique_elements(input_list):
    return dict(Counter(input_list))


def dict_to_sorted_string(element_counts):
    sorted_elements = sorted(element_counts.items(),
                             key=lambda x: x[1], reverse=True)
    result = ", ".join([f"{key} ({value})" for key, value in sorted_elements])
    return result


def create_binned_pcnts(percentages):
    bins = [0, 25, 50, 75, 90, 95, 100]

    percentages = [float(p) for p in percentages]

    labels = ['0%-25%', '25%-50%', '50%-75%', '75%-90%', '90%-95%', '95%-100%']
    if len(percentages) == 0:
        return ""

    binned = pd.cut(percentages, bins=bins, labels=labels,
                    right=True, include_lowest=True)
    counts = binned.value_counts().reindex(labels, fill_value=0)
    non_zero_counts = {label: count for label,
                       count in counts.items() if count > 0}
    result_str = ", ".join(
        [f"{label} ({count})" for label, count in non_zero_counts.items()])
    return result_str


def create_binned_seq_lens(numbers):
    # print("Numbers:", numbers)
    if len(numbers) == 0:
        return ""

    numbers = [int(n) for n in numbers]

    unique_counts = pd.Series(numbers).value_counts().to_dict()
    if len(unique_counts) < 6:
        return dict_to_sorted_string(unique_counts)

    bins = [0, 30, 100, 500, 1000, 3000, 5000, 10000, 1000000]
    labels = ['<30', '30-100', '100-500', '500-1000',
              '1000-3000', '3000-5000', '5000-10000', '>10000']

    binned = pd.cut(numbers, bins=bins, labels=labels,
                    right=True, include_lowest=True)
    counts = binned.value_counts().reindex(labels, fill_value=0)
    non_zero_counts = {label: count for label,
                       count in counts.items() if count > 0}
    total = sum(counts)
    result_str = ", ".join(
        [f"{label} ({count}, {count / total:.1%})" for label, count in non_zero_counts.items()])
    return result_str


def create_binnned_year(years):

    min_year = min(years)
    max_year = max(years)

    min_year = min_year // 10 * 10
    max_year = (max_year // 10 + (1 if (max_year % 10) else 0)) * 10 + 1

    bins = range(min_year, max_year, 10)
    labels = [f"{start}-{start+9}" for start in range(min_year, max_year - 1, 10)]

    binned = pd.cut(
        years,
        bins=bins,
        labels=labels,
        right=False,
        include_lowest=True)
    counts = binned.value_counts().reindex(labels, fill_value=0)
    non_zero_counts = {label: count for label,
                       count in counts.items() if count > 0}

    total_count = sum(non_zero_counts.values())
    percentages = {label: round((count / total_count) * 100, 1) for label, count in non_zero_counts.items()}

    result_str = ", ".join(
        [f"{label} ({count}) ({percentages[label]}%)" for label, count in non_zero_counts.items()])

    return result_str


def count_rev_sorter(value_count_list):
    return sorted(value_count_list, key=lambda x: (-x[-1], str(x[0])))


def alphabetical_sorter(value_count_list):
    return sorted(value_count_list, key=lambda x: x[0])


def int_sorter(value_count_list):
    return sorted(
        value_count_list,
        key=lambda x: int(x[0]) if (x[0] and x[0] != 'NA') else -1)


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

    return counter


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


def median_year(entry):

    year_list = str(entry).split(',')
    year_list = [i.replace('\u200b', '').strip().replace('–', '-') for i in year_list]

    def range_year(years):
        start, stop = years.split('-')
        return list(range(int(start), int(stop) + 1))

    year_list = [
        [year] if '-' not in year else range_year(year)
        for year in year_list
    ]

    year_list = [
        j
        for i in year_list
        for j in i
        if j
    ]

    year_list = [int(i) for i in year_list]

    return round(median(year_list)) if year_list else ''


def with_country(country):
    return 'Yes' if (country and country != 'NA') else 'No'


def dump_csv(file_path, table, headers=[], remain=True):

    file_path = Path(file_path)

    table_headers = []
    for rec in table:
        for key in rec.keys():
            if key not in table_headers:
                table_headers.append(key)

    if not headers:
        headers = table_headers
    else:
        remain_headers = [
            i
            for i in table_headers
            if i not in headers
        ]
        if remain:
            headers = headers + remain_headers
        table = [
            {
                k: v
                for k, v in i.items()
                if k in headers
            }
            for i in table
        ]

    file_path.parent.mkdir(exist_ok=True, parents=True)

    with open(file_path, 'w', encoding='utf-8-sig') as fd:
        writer = csv.DictWriter(fd, fieldnames=headers)
        writer.writeheader()
        writer.writerows(table)


def load_csv(file_path):
    records = []

    with open(file_path, encoding='utf-8-sig') as fd:
        for record in csv.DictReader(fd):
            records.append(record)
    return records


def format_counts_and_percentages(data_dict, total=None):

    if total is None:
        total = sum(data_dict.values())

    # Compute percentages with one decimal place
    percentages = {key: round((count / total) * 100, 1) for key, count in data_dict.items()}

    # Format output
    counts_formatted = "\n".join([f"{key}: {count}"
                    for key, count in count_rev_sorter(data_dict.items())])
    percentages_formatted = ", ".join([f"{key} ({percentage}%)"
                    for key, percentage in count_rev_sorter(percentages.items())])

    return counts_formatted, percentages_formatted


def yes_no(message, show_option=False, skip=False):

    if skip:
        return True

    option = input(f"{message} [y/n]:")

    option = option.lower()

    option = option == 'y'
    if show_option:
        if option:
            print('=> Choose Yes')
        else:
            print('=> Choose No')
    return option
