package config;

use strict;

use vars qw(%config);

our %config = (

	##### Tools - Rufus
	"bedtools" => "bedtools",

	##### Tools - Furious
	#"bedtools" => "/ifs/illumina/share/Utilities/softwareSL/BEDTools/bedtools-2.17.0/bin/bedtools",

	#"mm9" => "/ifs/illumina/share/Genomes/Mus_musculus/mm9/Fasta/Unmasked/mm9.fa.fai",  

	#"hg19" =>  "/ifs/illumina/share/Genomes/Homo_sapiens/hg19/Fasta/Unmasked/hg19.fa.fai",

	##### Genome data
	"mm9" => {
		"genome" => "/ifs/illumina/share/Genomes/Mus_musculus/mm9",
		"fa" => "/ifs/illumina/share/Genomes/Mus_musculus/mm9/Fasta/Unmasked/mm9.fa",
		"fai" => "/ifs/illumina/share/Genomes/Mus_musculus/mm9/Fasta/Unmasked/mm9.fa.fai",
	},

	"hg19" => {
		"genome" => "/ifs/illumina/share/Genomes/Homo_sapiens/hg19",
		"fa" => "/ifs/illumina/share/Genomes/Homo_sapiens/hg19/Fasta/Unmasked/hg19.fa",
		"fai" => "/ifs/illumina/share/Genomes/Homo_sapiens/hg19/Fasta/Unmasked/hg19.fa.fai",
	},

);

1;
