
def get_most_recent(years):
    if pd.isnull(years) or years == 'nan':  # Check for NaN or 'nan'
        return 'nan'
    # Split the string into individual years, clean them, and convert to integers
    year_list = [int(y.strip()) for y in str(years).split(',')]
    # Return the most recent year
    return max(year_list)


def remove_parenthesis(entry):
    # Remove everything in parentheses
    entry = re.sub(r"\s*\([^)]*\)", "", entry)
    # Extract the main part of the entry (excluding numbers)
    match = re.match(r"^[^0-9]+", entry)  # Correctly define 'match'
    # Return the matched portion if it exists; otherwise, return the cleaned entry
    return match.group().strip() if match else entry


def calculate_average_length(entry):
    # Split by commas and process each component
    components = entry.split(',')
    total_length = 0
    total_frequency = 0

    for component in components:
        # Remove leading/trailing whitespace
        component = component.strip()

        # Check for range
        range_match = re.match(r"(\d+)-(\d+)\s*\((\d+)\)", component)
        if range_match:
            start, end, freq = map(int, range_match.groups())
            avg_length = (start + end) / 2  # Average the range
            frequency = freq
        else:
            # Match single value with frequency
            value_match = re.match(r"(\d+)\s*\((\d+)\)", component)
            if value_match:
                value, freq = map(int, value_match.groups())
                avg_length = value
                frequency = freq
            else:
                # Skip invalid entries
                continue

        # Update totals
        total_length += avg_length * frequency
        total_frequency += frequency

    # Calculate weighted average

    return round(total_length / total_frequency) if total_frequency > 0 else 0


def most_frequent_range(entry):
    # Split by commas to process each range and frequency
    components = entry.split(',')
    max_frequency = 0
    most_frequent = None

    for component in components:
        # Match ranges and frequencies
        match = re.match(r"([\d%\-]+)\s*\((\d+)\)", component.strip())
        if match:
            range_, freq = match.groups()
            freq = int(freq)
            # Update the most frequent range
            if freq > max_frequency:
                max_frequency = freq
                most_frequent = range_

    return most_frequent

