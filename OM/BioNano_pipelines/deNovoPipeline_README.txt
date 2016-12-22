REAME for BioNano de novo assembly Pipeline

The de novo assembly pipeline is a set of python scripts which calls BioNano's RefAligner and Assembler to perform all the stages of de novo assembly and structural variation detection. BioNano's RefAligner and Assemblers are binary files used for alignment/refinement/merge and assembly, respectively.
The de novo assembly pipeline includes the following stages: 
	-Pairwise: pair-wise comparisons of all molecules, 
	-Assembly: construction of optical maps out of the pair-wise molecule alignments
	-Refinement: addition or removal of nick sites and re-position of nick sites based on molecules
	-Extension: extension of optical maps length by aligning molecules to the ends of optical maps
	-Merge: merging of the extended optical maps together if they overlap with each other
	-Final Refinement: refinment of the merged optical maps.
The structural variation includes the alignment of the assembled optical maps to a reference, and basing on the resulting alignments, the identification of insertion, deletion, inversion and translocation between the assembly and the reference.
The de novo assembly pipeline is run using the script pipelineCL.py