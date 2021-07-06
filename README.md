# Dinucleotide frequency
## Goal

This tools was used in the paper "MeCP2 is a microsatellite binding protein that protects CA repeats from nucleosome invasion", published in Science in 2021 (DOI: 10.1126/science.abd5581). 

CpA/CpG heatmaps were generated using the scripts available in this repository. The scripts generate matrices of dinucleotide counts per bin in regions of interest. TreeView v3.0 can then be used to plot heatmaps from generated matrices.

The tool requires BEDtools to run.

## Usage
Usage: countMotif.pl -w DIRNAME -b FILENAME -g STRING -m REGEXP -B INT

  -w DIRNAME or --workingDir=DIRNAME	 - Working directory (mandatory)
  -b FILENAME or --bed=FILENAME		 - Input bed file (mandatory)
  -m REGEXP or --motif=REGEXP 		 - regular expression standing for the motif to search for. eg: [ACTG] (mandatory)
  -g STRING or --genome=STRING 		 - Genome (mm9 or hg19) - (mandatory)
  -d INT or --distance=INT 		 - Distance arount TSS (default: 1500)
  --bin=INT 			   	 - bin size (default: 20)
