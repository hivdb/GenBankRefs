from dotenv import load_dotenv
from pathlib import Path
from openai import OpenAI, OpenAIError
import os
import time
import pandas as pd


load_dotenv()

TEMPLATE = """
    find the publication for this paper, and just give me the result:
    Title: {Title}
    Authors: {Authors}
    Journal: {Journal}
    Year: {Year}
    Accessions: {Accession}
"""


def build_prompts(system_prompt, metadata):

    prompt = TEMPLATE.format(**metadata)

    prompt = [
        {
            "role": "system",
            "content": system_prompt},
        {
            "role": "user",
            "content": prompt,
        }
    ]

    return prompt


def chat_openai(messages, model='gpt-4o', temperature=0):

    client = OpenAI(api_key=os.getenv('OPENAI_API_KEY'))

    try:
        response = client.responses.create(
            model=model,
            input=messages,
            tools=[{"type": "web_search_preview"}],
            temperature=temperature,
        )
    except OpenAIError as e:
        print(e)
        return []

    resp = [
        i
        for i in response.output
        if 'content' in i.model_fields
    ]
    return resp


def ask_ai(prompt):
    time.sleep(1)
    response = chat_openai(prompt)

    if len(response) == 0:
        return 'No response'
    else:
        return response[0].content[0].text


def using_ai_match(virus, genbank_unmatched, file_suffix, overwrite=False):
    virus_name = virus.full_name if virus.full_name else virus.name

    system_prompt = open(
        Path(__file__).resolve().parent / 'AI_match_template.txt').read().format(
            virus_name=virus_name)

    cache_file = virus.output_excel_dir / f'{virus.name}_{file_suffix}.xlsx'
    answers = []
    if cache_file.exists() and not overwrite:
        answers = pd.read_excel(cache_file)
        answer_map = {
            int(i['RefID']): i['AI_answer']
            for _, i in answers.iterrows()
        }
        for idx, row in genbank_unmatched.iterrows():
            genbank_unmatched.at[idx, 'AI_answer'] = answer_map.get(row['RefID'], '')

    for idx, row in genbank_unmatched.iterrows():
        if 'AI_answer' in row:
            continue
        # if 'Direct Submission' in row['Title']:
        #     continue
        # if 'Patent' in row['Journal']:
        #     continue

        prompt = build_prompts(system_prompt, {
            'Title': row['Title'],
            'Authors': row['Authors'],
            'Journal': row['Journal'] if row['Journal'].lower() != 'unpublished' else '',
            'Year': row['Year'] if row['Year'] else '',
            'Accession': row['accession']
        })
        # print(prompt)
        # raise

        answer = ask_ai(prompt)
        genbank_unmatched.loc[idx, 'AI_answer'] = answer

        answers.append({
            'RefID': row['RefID'],
            'Title': row['Title'],
            'Authors': row['Authors'],
            'Journal': row['Journal'] if row['Journal'].lower() != 'unpublished' else '',
            'Year': row['Year'] if row['Year'] else '',
            'Accession': row['accession'],
            'AI_answer': answer,
        })
        print(idx, 'done')

    pd.DataFrame(answers).to_excel(cache_file, index=False)

    return genbank_unmatched
