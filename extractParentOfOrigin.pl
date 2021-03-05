#!/usr/bin/perl

use warnings;
use strict;

###Script credit: Laurent Francioli (adapted from code provided at https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/2015-02-17-2014-09-09/5045-Step-order-for-PhaseBytransmission-and-ReadBackedPhasing)

my $usage = "<in.vcf> <in.dnms.intervals> child_id <out>\n";

#Check user input
die $usage if($#ARGV<3);

#Open files
open(IN,"<",$ARGV[0]) or die "ERROR: Could not open inupt VCF file: $ARGV[0]\n$usage";
open(DNM,"<",$ARGV[1]) or die "ERRPR: Could not open input DNM interval file: $ARGV[1]\n$usage";
open(OUT,">",$ARGV[3]) or die "ERROR: Could not open output file $ARGV[3]\n$usage";
my $child_id = $ARGV[2];

#Load DNM list
my %DNM;
while(<DNM>){
	chomp($_);
	$DNM{$_}=1;
}
close(DNM);

#Write output header
print OUT join("\t",("CHROM","POS","REF","ALT","CHILD_ID","HAP_ID","PAT_HAPLOTYPE_COUNT","MAT_HAPLOTYPE_COUNT","ORIGIN","LOCUS"))."\n";

#Hash to store header fields
my %HD;

#Hashes for storing results
#DNM IDs
my @dnms;
#DNM site-level info
#my %sites;
#DNM -> haplotype info
my %dnm_haps;
#Haplotype -> counts info
my %matCounts;
my %patCounts;
my %matLocus;
my %patLocus;

#Temporary vars for parsing
my @fields;
my @haps;
my $dnm_hap;

while(<IN>){
	#Skip description
	next if /^##/;
	chomp($_);
	@fields = split(/\t/,$_);

	#Parse header
	if($fields[0] eq "#CHROM"){
		die "Could not find child $child_id in input VCF file.\n$usage" if($_ !~ /$child_id/);
		for(my $i=0;$i<=$#fields;$i++){
			$HD{$fields[$i]}=$i;
		}
	}
	#Sites
	else{
		my $PhaseInconsistent = ($fields[$HD{"INFO"}] =~ /PhasingInconsistent/) ? 1 : 0;
		#Load child genotype fields in a hash
		my %child;
		@child{split(/:/,$fields[$HD{"FORMAT"}])} =  split(/:/,$fields[$HD{$child_id}]);
		my $gt = defined($child{"GT"}) ? $child{"GT"} : "";
		my $hp = (!$PhaseInconsistent && defined($child{"HP"}) && $child{"HP"} ne "." ) ? $child{"HP"} : "";

		#If haplotype information is available, use it 
		if($hp ne ""){
			@haps = split(/,/,$hp);	
			#If site is a DNM, add the haplotype of interest to the list
			if(defined($DNM{$fields[0].":".$fields[1]."-".$fields[1]})){
				#Find mutation haplotype
				$gt =~ /([01])\/([01])/;
				$dnm_hap = $1 == 1 ? $haps[0] : $haps[1];
				#Save site info in array (keep sorted) and use as key to find haplotype
				my $id = join("\t",(@fields[0,1,3,4],$child_id,$dnm_hap));
				push(@dnms,$id);
				$dnm_haps{$id} = $dnm_hap;

				#If mutation haplotype is not yet in the counts, initialize at 0
				$matCounts{$dnm_hap} = 0 if(!defined($matCounts{$dnm_hap}));
				$patCounts{$dnm_hap} = 0 if(!defined($patCounts{$dnm_hap}));
			}

			#If not but site is phased, count sites spanning the haplotype
			elsif($gt=~/\|/){
				$patCounts{$haps[1]}++;
				$patLocus{$haps[1]}=$fields[1];				
				$matCounts{$haps[0]}++;
				$matLocus{$haps[0]}=$fields[1];				
			}
		}
		#If the site is a DNM but the haplotype is unset or the phase inconsistent, report the site as unknown
		#Use a different coding for the site so that it cannot interfer with phase haplotypes
		elsif(defined($DNM{$fields[0].":".$fields[1]."-".$fields[1]})){
			my $id=join("\t",(@fields[0,1,3,4],$child_id,"NA"));
			push(@dnms,$id);
			$dnm_haps{$id}=$id;
			#Results are coded -- not the cleanest way:
			#-1=PhaseInconsistent
			#-2=Child homozygous
			#-3=Missing genotype
			#-4=Other
			if($PhaseInconsistent){
				$patCounts{$id}=-1;
				$matCounts{$id}=-1;
			}
			elsif($gt =~ /(0[\/\|]0)|(1[\/\|]1)/){
				$patCounts{$id}=-2;
				$matCounts{$id}=-2;
			}elsif($gt eq "."){
				$patCounts{$id}=-3;
				$matCounts{$id}=-3;
			}else{
				$patCounts{$id}=-4;
				$matCounts{$id}=-4;
			}
		}
	}
}

close(IN);

#Process all DNM sites and write results
foreach my $dnm(@dnms){
	print OUT join("\t",($dnm,$patCounts{$dnm_haps{$dnm}},$matCounts{$dnm_haps{$dnm}}));
	if($patCounts{$dnm_haps{$dnm}} > 0 && $matCounts{$dnm_haps{$dnm}} == 0){ print OUT "\tpaternal\t$patLocus{$dnm_haps{$dnm}}";}
	elsif($patCounts{$dnm_haps{$dnm}} == 0 && $matCounts{$dnm_haps{$dnm}} > 0){ print OUT "\tmaternal\t$matLocus{$dnm_haps{$dnm}}";}
	elsif($patCounts{$dnm_haps{$dnm}}==-1){print OUT "\tPhaseInconsistent";}
	elsif($patCounts{$dnm_haps{$dnm}}==-2){print OUT "\tChildHom";}
	elsif($patCounts{$dnm_haps{$dnm}}==-3){print OUT "\tMissingChildGT";}
	elsif($patCounts{$dnm_haps{$dnm}}==-4){print OUT "\tOtherPhaseProb";}
	else{ print OUT "\tunknown";}
	print OUT "\n";
}

close(OUT);