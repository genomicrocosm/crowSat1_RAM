README for BioNano hybrid scaffold Pipeline

The hybrid scaffold pipeline is a set of perl scripts which calls BioNano's RefAligner to perform all the stages of scaffolding. BioNano's RefAligner is a binary file used for alignment/refinement/merge.
The hybrid scaffold pipeline includes the following steps:
	-Conversion: conversion of the sequence FASTA file into CMAP format
	-Conflict identification and resolution: conflicts between the sequence and BioNano optical map assemblies can be identified and resolved automatically
	-Pairmerge: Merging of conflict-resolved sequences and optical maps to construct hybrid scaffolds
	-Final alignment: alignment of sequences and optical maps to the hybrid scaffolds
	-Conversion: conversion of the hybrid scaffolds from CMAP to AGP and FASTA formats
The hybrid scaffold pipeline is run using the hybridScaffold.pl script.