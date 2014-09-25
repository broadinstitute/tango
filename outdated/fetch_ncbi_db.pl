# This program automatically downloads a specified search entry from a specific database of NCBI 

#!/usr/bin/perl
use strict;
use Getopt::Long;
use LWP::Simple;

my %option = (
	h         => '',
	silent 		=> '',
	# I/O setting 
	db				=> '', 
	query			=> '',
	fetchdb		=> '',
	odir			=> '',
);

GetOptions(
	"h"							=> \$option{h},
	"silent"				=> \$option{silent},
	"db=s"					=> \$option{db},
	"query=s" 			=> \$option{query},
	"fetchdb=s"			=> \$option{fetchdb},
	"odir=s"				=> \$option{odir}
) || printHelp (); 

if ($option{h}) { printHelp();}

my $db = $option{db};
my $fetchdb = $option{fetchdb};
my @queries = split (',', $option{query});
my $odir = $option{odir};

unless ($odir && $option{query} && $db && $fetchdb) {
	print "\nErr: -db, -query, -fetchdb and -odir should be specified\n";
	printHelp ();
}

# E-Utilities base
my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

# get time of creation
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
++ $mon; #starting from 0
my $print_date = "-".$mday."-".$mon."-".$year;
                                             
foreach my $query (@queries) { # for each query of interest, query is in format of "viruses[organism]"

	my @fields = split ('\[', $query);
	
	# check of attempted download is human genome
	my $is_human_genome = 0;
	if (($fields[0] =~ /human/i) && ($db =~ /genome/i)){	
		unless ($option{silent}) {			
			print "\n\tHuman genome is being retrieved, only primary assembly will be included\n";
		}
		$is_human_genome = 1;	
	}
	

	# create output fasta file 
	my $output = $odir."/$db";
	$output .= "-".$fields[0].$print_date.".fasta";
	open (OUTPUT, ">$output") or die "cannot write to $output\n";
		
	$query .= "&retmax=10000000";
	
	# ----------- esearch ---------------
	my $url = $base . "esearch.fcgi?db=$db&term=$query";

	print "\n".$url."\n\n";   # debug
#		print get($url)."\n"; # debug
#		exit;
	
	my $num_sEntry;
	my $scounter = 0;
	my $num_genomes = 0;
	my @sEntries = split ("\n", get($url));
	
	foreach my $sEntry (@sEntries) {
	
		++ $scounter;
		
		if ($sEntry =~ /\<Count\>(\d+)\</) {
			$num_sEntry = $1;
			unless ($option{silent}) { print "\n$num_sEntry genomeID records writing to file $output\n"; }
		}	

		if ($scounter % 500 == 0) {
			unless ($option{silent}) { print "\t".$scounter." processed\n"; }
		}
						
		if ($sEntry =~ /\<Id\>(\d+)\</) { 
		
			my $genomeId = $1;				

			if ($db eq $fetchdb) {
					#-------- efetch -----------
					my $record = "";
					fetch_fasta (\$record, $base, $db, $genomeId, $is_human_genome);
					if ($record ne "") {
							print OUTPUT $record."\n";
							++ $num_genomes;							
					}								
			} else {
			
				#-------- elink -----------
				$url = $base . "elink.fcgi?dbfrom=$db&db=$fetchdb&id=$genomeId";   
			
				# completed chromosomes or other genomic sequences or plasmids and refseq
	#			my $lterm= "term=(gene+in+chromosome[prop]+OR+gene+in+genomic[prop]+OR+gene+in+plasmid[prop])+AND+srcdb+refseq[prop]";				

	#			my $lterm = "term=srcdb+ddbj/embl/genbank[prop]";  # INSDC
			
	#			$url .= "&$lterm";

	#				print $url."\n";   # debug
	#				print get($url)."\n"; # debug
	# 			exit;

				my @lEntries = split ("\n", get($url));
			
				my $start_parsing = 0; # LinkSetDb tag is found 
			
				foreach my $lEntry (@lEntries) {
			
					if ($lEntry =~ /LinkSetDb/) { $start_parsing = 1; } # skip Id that corresponding to genomeID
				
					if (($lEntry =~ /\<Id\>(\d+)\</) && ($start_parsing == 1)) {
				
						my $nucId = $1;

						#-------- efetch -----------
						my $record = "";
						fetch_fasta (\$record, $base, $fetchdb, $nucId, $is_human_genome);
						if ($record ne "") {
								print OUTPUT $record."\n";
								++ $num_genomes;							
						}
						
#						$url = $base . "efetch.fcgi?db=$fetchdb&id=$nucId";
#						$url .= "&rettype=fasta&retmode=text";
#						my $record = get ($url); 
#						$record =~ /(.+?)\n(.+)/;
#						my $header = $1;
#						my $seq = $2;
#						if ($seq ne "") {   # do not print empty record
#							while ($record =~ /\n$/) { chomp ($record); } # remove all \n from end of the record										
#							if ($is_human_genome == 1) {
#								if ($header =~ /Primary/) {								
#									print OUTPUT $record."\n";														
#									++ $num_genomes;
#								}
#							} else {
#								print OUTPUT $record."\n";
#								++ $num_genomes;
#							}
#						} # if ($seq ne "") 
						
					} #	if (($lEntry =~ /\<Id\>(\d+)\</) && ($start_parsing == 1)) 
				
				} # foreach my $lEntry (@lEntries) 
				
			} # if ($db ne $fetchdb) 						
			
		} # 			if ($sEntry =~ /\<Id\>(\d+)\</)  

	}	#	foreach my $sEntry (@sEntries) {

	close (OUTPUT);	

	unless ($option{silent}) { print "\n\t".$num_genomes." sequencesv written\n"; }
			
} # foreach my $query (@queries)
	
### include this code for ESearch-EFetch

#post the efetch URL

################################################################################################################################
# sub-functions 
#################################################################################################################################
sub printHelp {
		print "\n----------------------------------------------------------------------------\n";
		print "Brief: Retrieve sequences from NCBI database\n";
		print "\nusage: ./$0 -db [searchdb] -query [q1,q2...] -fetchdb [fetchdb] -odir [odir]\n\n";
		print "-h: help; for cmd options\n";
		print "-silent: no screen output for programs called\n\n";
		print "I/O setting\n";
		print "-db: comma separated NCBI databases\n";
		print "-query: comma separated esearch terms, e.g. viruses[organism]\n";
		print "-fetchdb: database to fetch sequence from \n";
		print "-odir: output DIR; each pair of db+term will become one output file\n";
		print "\nExample: perl fetch_ncbi_db.pl -db genome -query viruses\[organism\],human\[organism\],fungi\[organism\],bacteria\[organism\] 
		\n\t-fetchdb nuccore -odir /seq/viral/analysis/xyang/FUO/curated_database/\n";
		print "\n----------------------------------------------------------------------------\n";
		exit;
}


sub fetch_fasta {
	my $ref_record = shift;
	my $base = shift;
	my $db = shift;
	my $id = shift;
	my $is_human_genome = shift;
	
	#-------- efetch -----------
	my $url = $base . "efetch.fcgi?db=$db&id=$id";
	$url .= "&rettype=fasta&retmode=text";

	$$ref_record = get ($url); 
	$$ref_record =~ /(.+?)\n(.+)/;
	my $header = $1;
	my $seq = $2;

	if ($seq ne "") {   # do not print empty record

		while ($$ref_record =~ /\n$/) { chomp ($$ref_record); } # remove all \n from end of the record
						
		if ($is_human_genome == 1) {
			unless ($header =~ /Primary/) {			
				$$ref_record = "";
			}
		} 
	} else {
		$$ref_record = "";
	}
} # fetch_fasta
