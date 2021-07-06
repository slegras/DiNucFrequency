package countMotif;

use strict;
use warnings;
use Exporter;
use POSIX;

our @ISA= qw(Exporter);
our @EXPORT= qw(&extendWindow &binning &getFastafromBED &countMotif &formatIn &formatWindow);

use config;
use vars qw(%config );
*config=\%config::config;

sub binning{
	my ($bed, $bin, $motif, $out, $log) = @_;

	## Build up the command line
	## Sur cette commande, je reviens toujours un peu en avant de la fin pour compter les
	## motifs qui sont à cheval des fenêtres
	my $step = $bin-(length($motif)-1);
	my $cmdline = "$config{'bedtools'} makewindows -b $bed -w $bin -s $step -i srcwinnum > $out";

	my $title = "Binning regions.";

	&runCmdLine($cmdline, $log, $title);

	return 1;
}

sub countMotif{
	my ($fa, $distance, $bin, $motif, $out, $log) = @_;

	if ($fa =~ /.gz$/){
		`gunzip $fa`;
		$fa =~ s/\.gz$//;
	}

	open(FA, "<".$fa) or die "Cannot open fasta file $fa: $!";
	open(OUT, ">".$out) or die "Cannot create output file: $!";

	## printing headers
	print OUT "Gene";
	my $i=0;
	for($i = -$distance; $i<= $distance-$bin-(length($motif)-1); $i+=$bin-(length($motif)-1)){
		print OUT "\t".$i;
	}
	print OUT "\t".$i if($i+$bin < $distance);
	print OUT "\t".$distance."\n";


	## Writing the number of CpG per sequences
	my $seq;
	my $nb_motif=0;
	my $id;
	while (<FA>){
		## On vire le \n de fin de ligne
		chomp;

		## Lorsque l'on lit une nouvelle entête de séquence
		if(/^>/){

			#print "$_\n";

			s/^>//;
			s/_[0-9]+$//;
			s/__$//;
			s/_\+$//;

			## Si on a déjà lu une séquence aupavant
			## On compte le nombre de motif et on print ce nombre
			if(defined $id){
				## on regarde dans la séquence et on compte
				$nb_motif++ while($seq =~ /$motif/g);

				## On détermine la séquence reverse complémentaire du motif
				my $revComp = &revComp($motif);
				## On compte le nombre d'occurence pour le revcomp du motif
				$nb_motif++ while($seq =~ /$revComp/g);

				print OUT "\t".$nb_motif;

				# print "$id : $_\n";

				unless(/^$id$/){
					# print "No\n";
					print OUT "\n";
					print OUT $_;
				}

				## On réinitialise la sequence et le nombre de motif
				## car on commence une nouvelle sequence
				$nb_motif = 0;
				$seq = "";

			}else{
				#print "troupe\n";
				print OUT $_;
			}

			## On enregistre le nom de la séquence en cours de lecture
			$id = $_;

		}
		## Lorsque l'on lit une séquence nucléotidique, on concatène la séquence aux
		## autre ligne de séquence vu précédemment (pour la même entête)
		else{
			$seq .= $_;
		}

	}
	## Pour la dernière ligne avec la séquence
	$nb_motif = 0;
	## on regarde dans la séquence et on compte
	$nb_motif++ while($seq =~ /$motif/g);

	## On détermine la séquence reverse complémentaire du motif
	my $revComp = &revComp($motif);
	## On compte le nombre d'occurence pour le revcomp du motif
	$nb_motif++ while($seq =~ /$revComp/g);
	print OUT "\t".$nb_motif."\n";

	close FA;
	close OUT;

	`gzip $fa`;

	return 1;
}

sub extendWindow{
	my ($bed, $distance, $genome, $out, $log) = @_;

	## Build up the command line
	my $cmdline = "$config{'bedtools'} flank -g $config{$genome}{fai} -b $distance -i $bed > $out";

	my $title = "Extending input coordinates of $distance each side.";

	&runCmdLine($cmdline, $log, $title);

	return 1;
}

sub formatIn{
	my ($in, $out, $log) = @_;

	open(IN, "<".$in) or die "Cannot open input file $in: $!";
	open(OUT, ">".$out) or die "Cannot open output file $out: $!";

	while(<IN>){
		#s/\r\n//;
		chomp;
		my @tab = split "\t";
		my $center = $tab[1]+($tab[2]-$tab[1])/2;
		$center = floor($center);
		print OUT $tab[0]."\t";
		print OUT $center."\t";
		print OUT $center;
		print OUT "\t";
		#print OUT $tab[3]."_".$tab[4]."_".$tab[5]."\t";
		# print OUT $tab[0]."_".$tab[1]."_".$tab[2]."_".$tab[4]."_".$tab[5]."\t";
		print OUT $tab[0]."_".$tab[1]."_".$tab[2]."_".$tab[5]."\t";
		print OUT $tab[4]."\t";
		print OUT $tab[5]."\n";
	}

	close IN;
	close OUT;

	return 1;
}

sub formatWindow{
	my ($in, $out, $log) = @_;
	open(IN, "<".$in) or die "Cannot open input file $in: $!";
	open(OUT, ">".$out) or die "Cannot open output file $out: $!";

	while(<IN>){
		chomp;
		my @tab = split "\t";

		my ($gene, $strand) = $tab[3] =~ /(.+)_([\+-])_[0-9]+$/;

		print OUT $tab[0]."\t";
		print OUT $tab[1]."\t";
		print OUT $tab[2];
		print OUT "\t";
		print OUT $tab[3]."\t";
		print OUT $gene."\t";
		print OUT $strand."\n";
	}

	close IN;
	close OUT;


	return 1;
}

sub getFastafromBED{
	my ($bed, $genome, $out, $log) = @_;

	## Build up the command line
	my $cmdline = "$config{'bedtools'} getfasta -fi $config{$genome}{fa} -bed $bed -name -fo - | \
		perl -ne 's/(::.*)//; s/[:-]/_/g; print uc(\$_)' | gzip > $out.gz";

	my $title = "Extracting nucleotide sequences";

	&runCmdLine($cmdline, $log, $title);

	return 1;

}

sub revComp() {
        my %baseComp = (
                'A' => 'T',
                'C' => 'G',
                'T' => 'A',
                'G' => 'C',
                '[' => ']',
                ']' => '['
        );

        my $seq     = $_[0];
        my @rev     = split "", reverse($seq);
        my $revComp = "";
        foreach my $base (@rev) {
                $revComp .= $baseComp{$base} if($baseComp{$base});
        }
        return ($revComp);
}

sub runCmdLine {
	my ($cmdline, $log, $title) = @_;

	## Run the command line and output into the log file
	print $log "\n##############\n" ;
	print $log "## Start analysis: ".`date`."\n";
	print $log "## $title. \n\n";
	print $log "$cmdline\n";
	print $log "#<--- Output: --------------------------------------------------\n";

	my $result = `$cmdline 2>&1`;
	print $log $result."\n";

	print $log "\n";

	print $log "#--------------------------------------------------------------->\n";
	print $log "## End of analysis: ".`date`."\n";

	return 1;
}

1;
