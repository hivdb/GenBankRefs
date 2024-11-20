import pandas as pd
import numpy as np
import re
import os
from collections import Counter


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
        processed_names.append(processed_name)

    processed_names = ', '.join(processed_names)
    return (processed_names)


def extract_year_from_journal(text):
    match = re.search(r'\((\d{4})\)', text)
    if match:
        return match.group(1)
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
        return "NA"


# Each list is a string in which items are separated by ', '
def combine_items_in_different_lists(lists):
    unique_dict = {}
    for string in lists:
        if isinstance(string, int):
            string = str(string)
        items = string.split(', ')
        for item in items:
            unique_dict[item] = 1
    list_unique_items = list(unique_dict.keys())
    list_unique_items = sorted(list_unique_items)
    unique_items = ', '.join(list_unique_items)
    if unique_items.startswith(","):
        unique_items = unique_items[1:]
    unique_items = unique_items.strip()
    return unique_items


# Convert {key: [values], key: [values]] to a list of sets comprising {key, values}
# In this program, each key is the index of a row and the values contain
# one or more indexes of rows that share some property
# The subset flag code ensures sets do not share items with one another
def convert_dict_to_list_of_sets(matched_indexes):
    list_of_sets = []

    for row_i, row_j_list in matched_indexes.items():
        matched_rows = set()
        matched_rows.add(row_i)
        matched_rows.update(set(row_j_list))

        find_same_group = False
        for merged_indexes in list_of_sets:
            for row in matched_rows:
                if row in merged_indexes:
                    find_same_group = True
                    merged_indexes.update(matched_rows)

        if not find_same_group:
            list_of_sets.append(matched_rows)

    list_of_rows = [
        j
        for i in list_of_sets
        for j in i
    ]

    return (list_of_sets, list_of_rows)


def get_pcnt_authors_overlap(authors1, authors2):
    set1 = set(authors1.split(', '))
    set2 = set(authors2.split(', '))
    shared_set = set1 & set2
    combined_set = set1 | set2
    pcnt = len(shared_set) / len(combined_set)
    return pcnt


def get_pcnt_shared_accessions(stringlist1, stringlist2):
    list1 = stringlist1.split(', ')
    list2 = stringlist2.split(', ')
    pcnt_shared_accessions = len(set(list1) & set(
        list2)) / len(set(list1) | set(list2))
    return pcnt_shared_accessions


def get_pcnt_shared_stems(list1, list2, stem_length):
    acc_num_stem_list1 = []
    acc_num_stem_list2 = []
    for acc_num in list1:
        acc_num_stem_list1.append(acc_num[:stem_length])
    for acc_num in list2:
        acc_num_stem_list2.append(acc_num[:stem_length])
    set_acc_num_stem1 = set(acc_num_stem_list1)
    set_acc_num_stem2 = set(acc_num_stem_list2)
    pcnt_shared_stems = len(set_acc_num_stem1 & set_acc_num_stem2) / \
        len((set_acc_num_stem1 | set_acc_num_stem2))
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
    result_str = ", ".join(
        [f"{label} ({count})" for label, count in non_zero_counts.items()])
    return result_str
