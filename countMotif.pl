#! /usr/bin/perl

BEGIN {
	use File::Basename;
	push( @INC, dirname(__FILE__) );
}

use strict;
use warnings;
use Getopt::Long;
use lib qw(.);
use countMotif;

use config;
use vars qw(%config );
*config=\%config::config;

## Date : 25 Sep 2014
## Author : Stephanie Le Gras

## Objectives :

my $bed;
my $workingDir;
my $distance;
my $genome;
my $bin;
my $motif;

my $num_arg  = scalar @ARGV;

my $result = GetOptions(
	"b|bed=s"      => \$bed,
	"w|workingDir=s" => \$workingDir,
	"d|distance=s"  => \$distance,
	"g|genome=s" => \$genome,
	"bin=s" => \$bin,
	"m|motif=s" => \$motif,

);

my $usage = <<END;

Usage: $0 -w DIRNAME -b FILENAME -g STRING -m REGEXP -B INT

  -w DIRNAME or --workingDir=DIRNAME	 - Working directory (mandatory)
  -b FILENAME or --bed=FILENAME		 - Input bed file (mandatory)
  -m REGEXP or --motif=REGEXP 		 - regular expression standing for the motif to search for. eg: [ACTG] (mandatory)
  -g STRING or --genome=STRING 		 - Genome (mm9 or hg19) - (mandatory)
  -d INT or --distance=INT 		 - Distance arount TSS (default: 1500)
  --bin=INT 			   	 - bin size (default: 20)

END

die $usage if (@ARGV);
die $usage if (! $result);
die $usage if ( $num_arg < 4 );

### Checking settings
die "No input file.\n" unless ($bed);
die "No working dir.\n" unless ($workingDir and -d $workingDir);
die "Genome should be either mm9 or h19.\n" unless($genome =~ /(hg19)|(mm9)/);
die "Distance should be an integer\n" if(defined $distance and $distance !~ /^[0-9]+$/);

$distance = 1500 unless(defined $distance);
$bin = 20 unless(defined $bin);

##################################################################
################################# Creating a tmp dir
my $tmp = "$workingDir/tmp";
if(!-d $tmp){
	mkdir $tmp;
}

##################################################################
################################# Creating a log file
## Creating the LOG directory if it doesn't exist
if(!-d "$workingDir/LOG"){
	mkdir "$workingDir/LOG";
}

## Getting the log file name
if ( -f "$workingDir/".basename($0).".log"){
	unlink "$workingDir/".basename($0).".log";
}

my $i = 1;
my $logName = "$workingDir/LOG/".basename($0).".$i.log";

while ( -f $logName ){
	$i++;
	$logName = "$workingDir/LOG/".basename($0).".$i.log";
}

symlink $logName, "$workingDir/".basename($0).".log";

open(my $log, ">".$logName) or die "Cannot create log file : $!";


##################################################################
################################# MAIN
print $log "################ Starting Analysis\n";
print $log "##".`date`;
print $log "\n";

######### formating input file
my $out = "$tmp/1-formated.bed";
if (! -f $out){
	&formatIn( $bed, $out, $log );
}

######### Limited the analysis to one copy of each gene
$out = "$tmp/2-uniq_formated.bed";
if (! -f $out){
	`sort -u "$tmp/1-formated.bed" > $out`
}

######### Extending input coordinates
$out = "$tmp/3-window.bed";
if (! -f $out){
	&extendWindow( "$tmp/2-uniq_formated.bed", $distance, $genome, $out, $log );
}

######### Binning regions
$out = "$tmp/4-window_binned.bed";
if (! -f $out){
	&binning("$tmp/3-window.bed", $bin, $motif, $out, $log);
}

######### Formating file for fasta extraction (need strand addition)
$out = "$tmp/5-formated_window_binned.bed";
if (! -f "$out"){
	&formatWindow("$tmp/4-window_binned.bed", $out, $log);
}

######### Extracting Fasta sequences
$out = "$tmp/6-window_binned.fa";
if (! -f "$out.gz"){
	&getFastafromBED("$tmp/5-formated_window_binned.bed", $genome, $out, $log);
}

######### Counting Motifs
$out = "$tmp/7-window_binned_$motif.tsv";
if (! -f $out){
	&countMotif("$tmp/6-window_binned.fa.gz", $distance, $bin, $motif, $out, $log);
}
