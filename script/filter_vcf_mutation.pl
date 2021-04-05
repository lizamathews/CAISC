#!/bin/perl

# this perl script filters vcf files and generates new vcf files with 
# mutations that only appear in at least a certain number of cells

# ask for input VCF file name
print "Enter the name of the VCF you want to filter: ";
my $inputvcf = <STDIN>;
chomp $inputvcf;

# ask for mutation filtering threshold
print "Enter the minimum number of cells that mutations should appear in: ";
my $mutations = <STDIN>;
chomp $mutations;

# create output VCF file name
my $outputvcf = substr $inputvcf, 0, -4;
my $outputvcf = ">${outputvcf}_${mutations}mutations.vcf";

# filter for mutations that only appear in at least a certain number of cells
open(INFILE, $inputvcf)||die"";
open(OUTFILE, $outputvcf)||die"";
for($i=0; $i<126;$i++){
	$line=<INFILE>;
	print OUTFILE $line;
}
while($line=<INFILE>){
	chop($line);
	@items=split("\t", $line);
	$mut=0;
	for($i=9; $i<@items;$i++){
		if(@items[$i] ne "./."){
			$mut++;
		}
	}
	if($mut>=$mutations){
		print OUTFILE $line . "\n";
	}
}
close(OUTFILE);
close(INFILE);