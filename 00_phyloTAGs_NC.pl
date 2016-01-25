#!/usr/bin/perl -w
use strict;
use Getopt::Long;

###########################################################################################
# Script Name: 00_phyloTAGs_NC.pl
# Script Version:1.0
# Author: A. Caro-Quintero
# Date: July 2015
# 
#
# Description: This script facilitates the identification of the NC numbers, it takes 
# as a parameter the genus and/or species that will be used as a reference and retrieves 
# all possible NC numbers associated to the specified parameters. It can also print the 
# organisms and identity of the taxonomic group that will be target for primer designing. 
# To do this, the 16S rDNA identity range has to be specified. This information  
# can be used to established, (1) if the required reference is in the database, (2) the 
# number of related organisms within the targeted taxonomic range and (3) their 
# divergence with respect to the reference. All this information can be used to establish
# the level of representation of the targeted taxonomic group and the posible biases on the
# primer design.
#
# Parameters:
#
# --genus
# --species
# --NC
# 
#
###########################################################################################
my$h=0;
my$genus="";
my$species="";
my$NC_file;
my$k="";

GetOptions("genus:s" => \$genus, "species:s" => \$species, "NC:s" => \$NC_file, "help"=>\$h);

if ($genus ne $k)
{
$k=1;
}


if ($h > 0)
{
    print "\n\n Script: 00_phyloTAGs_NC.pl \n Author: Alejandro Caro-Quintero \n Version: 1.0 \n\n Description:\n -genus (genus of the reference, required)\n -species (species of the reference, optional)\n -NC (file with NC numbers map to names)\n\n\n";


}
elsif ($k>0) 
{
    print "\n ############################";
    print "\n  List of identifiers";
    print "\n ############################\n";
    open (FILE,"$NC_file")||die ("Could not open file $NC_file \n");
    while(<FILE>)
    {
	chomp;
	my$line=$_;
	if($line=~ /$genus/)
	{
	    if($line=~ /$species/)
	    {
		
		print " $line\n";
		
	    }
     	}
    }

    print "\n\n Do you want to find the NC numbers and 16S rDNA identity values \n\n associated to the reference? (y/n):";
    my$Y_N=<STDIN>;
    if($Y_N=~ /y|Y|yes|YES/)
    {
	print "\n 1. Please specified which identifier you want to use : ";
	my$NC_id=<STDIN>;
	chop($NC_id);
	print "\n 2. Please specified the distance file ? :";
	my$distance=<STDIN>;
	print "\n 3. Lower identity range? :";
	my$low=<STDIN>;
	chop($low);
	print "\n 4. Higher identity range? :";
	my$high=<STDIN>;
	chop($high);
	print "\n 5. Output file :";
	my$output=<STDIN>;
	print "\n\n";
	open(OUT, ">>$output");		
	open(DIST, $distance)||die("Could not open file $distance\n");
	while(<DIST>)
	{
	    chomp;
	    my$line2=$_;
	    my@array=split("\t",$line2);
	    if($array[0]=~ /$NC_id/)
	    {
		if($array[2] >$low && $array[2]<$high)
		{
		    print OUT "$line2\n";	
		}		
	    }
	print "\n #################\n";	
	print " ## All done!\n";
	print " #################\n";    
	}
    }	
    elsif ($Y_N=~ /n|N|no|NO/)
    {
	print "\n #################\n";	
	print " ## All done!\n";
	print " #################\n";
    }
    else 
    {
	print "\n ##########################\n";		
	print " ## Wrong parameter used!\n";
	print " ##########################\n";
    }
}
