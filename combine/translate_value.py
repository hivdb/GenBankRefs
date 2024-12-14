from statistics import median


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
