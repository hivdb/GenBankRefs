import importlib.util
import importlib
from pathlib import Path


# This is called by the main function to allow users to select the virus
def select_virus():
    viruses_list = (
        ('Orthonairovirus haemorrhagiae', 'CCHF'),
        ('Henipavirus nipahense', 'Nipah'),
    )
    for index, (cname, name) in enumerate(viruses_list):
        print(f"{index + 1}.", cname, f"({name})")

    virus_id = input('Please select a virus by ID: ')
    assert virus_id.isdigit(), 'Virus not found'
    assert int(virus_id) <= len(viruses_list), 'Virus not found'

    return viruses_list[int(virus_id) - 1][-1]


def load_virus_obj(virus):
    """
        virus_conf folder contains specific configuration (virus config) for a virus.
        If this file doesn't exist, the default.py file will be used.

        A virus config contains the input and output file path, and functions
        for cleaning the data or functions for running blast.
    """
    virus_conf_path = f'viruses/{virus}.py'
    if not Path(virus_conf_path).exists():
        spec = importlib.util.spec_from_file_location(
            'viruses.default', 'viruses/default.py')
    else:
        spec = importlib.util.spec_from_file_location(
            f'viruses.{virus}', virus_conf_path)

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not Path(virus_conf_path).exists():
        module.set_virus(virus)

    return module
