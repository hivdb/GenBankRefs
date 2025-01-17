
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path


def dump_fasta(file_path, seq_list):

    seq_list = [
            SeqRecord(
                id=str(k),
                seq=Seq(v),
                description=''
            )
            for k, v in seq_list.items()
        ]

    file_path = Path(file_path)
    file_path.parent.mkdir(exist_ok=True, parents=True)

    SeqIO.write(
        seq_list,
        str(file_path),
        'fasta')


def load_fasta(file_path):
    seq_list = {}
    for i in SeqIO.parse(str(file_path), 'fasta'):
        seq_list[i.name] = str(i.seq)

    return seq_list
