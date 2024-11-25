from statistics import median


def categorize_host_specimen(df, host_col, specimen_col):

    for index, row in df.iterrows():
        host = row[host_col].lower()
        specimen = row[specimen_col].lower()

        updated_host = []
        updated_specimen = []

        if 'homo sapiens' in host:
            updated_host.append('Homo sapiens')

        for a in ['animal', 'sheep', 'cattle', 'goat', 'mouse', 'boar', 'hare',
                  'livestock', 'cow', 'sheep', 'camel', 'monkey', 'deer',
                  'rodent', 'serow']:
            if a in host:
                updated_host.append('animal')

        if 'tick' in host:
            updated_host.append('Ticks')

        if 'tick' in specimen:
            updated_host.append('Ticks')

        if 'NA'.lower() == host or not host:
            updated_host.append('NA')

        updated_host.append('Other')
        updated_host.append(host)

        for i in ['serum', 'blood', 'plasma', 'sera']:
            if i in specimen:
                updated_specimen.append('blood')

        for s in ['brain', 'spleen', 'nasal swab']:
            if s in specimen:
                updated_specimen.append(s)

        if not specimen or specimen == 'NA'.lower():
            updated_specimen.append('NA')

        updated_specimen.append('Other')
        updated_specimen.append(specimen)

        df.at[index, 'CleanedHost'] = updated_host[0]
        df.at[index, 'CleanedSpecimen'] = updated_specimen[0]

    return df


def median_year(entry):

    year_list = entry.split(',')
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


def translate_country(country):
    return 'Yes' if (country and country != 'NA') else 'No'


def translate_gene(gene):
    return gene if gene in ['L', 'S', 'M'] else ('NA' if not gene or gene == 'NA' else 'Other')


def translate_specimen(specimen):
    specimen = specimen.lower()
    for i in ['serum', 'blood', 'plasma', 'sera']:
        if i in specimen:
            return 'blood'

    for s in ['brain', 'spleen', 'nasal swab']:
        if s in specimen:
            return s

    if not specimen or specimen == 'NA'.lower():
        return 'NA'

    return 'Other'


def translate_hosts(host):
    # Genbank pre process already did some translation,
    # here only use the processed one host result, and aggregate to some categories
    host = host.lower()
    if 'homo sapiens' in host:
        return 'Homo sapiens'

    for a in [
            'animal', 'sheep', 'cattle', 'goat', 'mouse', 'boar', 'hare',
            'livestock', 'cow', 'sheep', 'camel', 'monkey', 'deer',
            'rodent', 'serow']:
        if a in host:
            return 'animal'

    if 'tick' in host:
        return 'Ticks'

    if 'NA'.lower() == host or not host:
        return 'NA'

    return 'Other'
