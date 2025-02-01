from dotenv import load_dotenv
from pathlib import Path
from openai import OpenAI, OpenAIError
import os
import time
from Utilities import load_csv
from Utilities import dump_csv


load_dotenv()

TEMPLATE = """
    can you help me find the publication for this paper:
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
        response = client.chat.completions.create(
            model=model,
            messages=messages,
            temperature=temperature,
        )
    except OpenAIError as e:
        print(e)
        return []

    return response.choices


def ask_ai(prompt):
    time.sleep(1)
    response = chat_openai(prompt)

    if len(response) == 0:
        return 'No response'
    else:
        return response[0].message.content


def using_ai_match(virus, genbank_unmatched):
    system_prompt = f"""
        You are an AI expert in virology with specialized
        knowledge in {virus.name} virus, including its transmission,
        pathogenesis, and clinical management.
        You are also skilled in systematic review methods,
        capable of synthesizing research findings
        and providing evidence-based insights.
    """

    cache_file = virus.output_excel_dir / 'AI_cache.csv'
    answers = []
    if cache_file.exists():
        answers = load_csv(cache_file)
        answer_map = {
            int(i['RefID']): i['AI_answer']
            for i in answers
        }
        for idx, row in genbank_unmatched.iterrows():
            genbank_unmatched.at[idx, 'AI_answer'] = answer_map.get(row['RefID'], '')

    for idx, row in genbank_unmatched.iterrows():
        if 'AI_answer' in row:
            continue
        if 'Direct Submission' in row['Title']:
            continue
        if 'Patent' in row['Journal']:
            continue

        prompt = build_prompts(system_prompt, {
            'Title': row['Title'],
            'Authors': row['Authors'],
            'Journal': row['Journal'],
            'Year': row['Year'],
            'Accession': row['accession']
        })
        # print(prompt)

        answer = ask_ai(prompt)
        genbank_unmatched.at[idx, 'AI_answer'] = answer

        answers.append({
            'RefID': row['RefID'],
            'AI_answer': answer,
        })

    dump_csv(cache_file, answers)

    return genbank_unmatched
