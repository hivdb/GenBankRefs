#!/bin/bash

rm -rf output_L
ViralMSA.py -e hivdbteam@list.stanford.edu -s L_isolates.fasta -o output_L -r L_ref_na.fasta --omit_ref
iqtree2 -s L_isolates.fasta.aln -m GTR+G4+F -bb 1000 -nt AUTO

rm -rr output_N
ViralMSA.py -e hivdbteam@list.stanford.edu  -s N_isolates.fasta -o output_N -r N_ref_na.fasta --omit_ref
iqtree2 -s N_isolates.fasta.aln -m GTR+G4+F -bb 1000 -nt AUTO
