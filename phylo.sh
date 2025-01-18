#!/bin/bash

rm -f output_L
ViralMSA.py -s L_isolates.fasta -o output_L -r L_ref_na.fasta
iqtree2 -s L_isolates.fasta.aln -m GTR+G4+F -bb 1000 -nt AUTO

rm -r output_N
ViralMSA.py -s N_isolates.fasta -o output_N -r N_ref_na.fasta
iqtree2 -s N_isolates.fasta.aln -m GTR+G4+F -bb 1000 -nt AUTO
