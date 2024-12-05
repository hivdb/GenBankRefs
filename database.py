import pandas as pd
import sqlite3
from functools import reduce


def dump_table(db_file, table_name, table):
    conn = sqlite3.connect(str(db_file))

    table.to_sql(table_name, conn, if_exists='replace', index=False)


def load_table(db_file, table_name):
    conn = sqlite3.connect(str(db_file))
    return pd.read_sql(f'SELECT * FROM {table_name};', conn)


def load_tables(db_file, table_name_list):
    tables = []
    for i in table_name_list:
        tables.append(load_table(db_file, i))

    return tables


def split_table(table_config, dataframe):
    tables = []
    for table_name, columns in table_config:
        table = dataframe[columns]
        tables.append((table_name, table))

    return tables


def merge_table(tables, key):
    return reduce(lambda x, y: pd.merge(x, y, on=key, how='inner'), tables)
