#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#**********************************************************************************************************
# Script Name: 01_phyloTAGs_align.pl
# Script Version: 1.0
# Author: A.Caro-Quintero
# Date: August 2015
#
# Description: This script extracts and aligns by codons gene sequences, e.g., gyrB genes, related to a defined 
# reference with an specific divergence range. Evolutionary relatedness is provided in a separate file, 
# were pairwise 16S rDNA identity of all available completed sequences had been calculated, other distance 
# measurements such as the Average Aminoacid Identity (AAI) can also be used as as long as the files have the 
# same format. If the extraction of related genes is not required the user can specify the "no_distance" option, 
# in this case only the alignment by codons will be done. 
#
# The alignments by codons is done with ClustalW (Larkin et al. 2007) and the PAL2NAL script (Suyana et al. 2006). 
# 
# Parameters: 
#
# --in
# --no_distance
# --distance
# --reference
# --idlow
# --idhigh
#
# References:
#
# Mikita Suyama, David Torrents, and Peer Bork. PAL2NAL: robust conversion of protein 
# sequence alignments into the corresponding codon alignments.(2006) Nucleic Acids Res. 34, 
# W609-W612.
#
# Larkin, Mark A., et al. Clustal W and Clustal X version 2.0.(2007) Bioinformatics 23.21: 2947-2948.
#
#############################################################################################################

my$path=`pwd`;
chop($path);
my$alignment; #one line tab alignment with out the ">" sign having NC_# as identifiers
my$SrRNA_id;
my$reference;
my$lowest=0;
my$highest=0;
my$h=0;
my$pass=0;
my$nodist=0;

# other variables

my%files;
my%positions;
my$clustalw="clustalw";# or path to ClustalW (e.g ~/ClustalW/clustalw1.83/clustalw)
my$pal2nal="perl pal2nal\.v14/pal2nal\.pl";#~/pal2nal/pal2nal\.pl

# if not specified the values will be taken by default 
GetOptions("in:s" => \$alignment, "no_distance" => \$nodist, "distance:s" => \$SrRNA_id, "reference:s" => \$reference, "idlow:s"=>\$lowest, "idhigh:s"=>\$highest, "help"=>\$h );

# printing help 
if ($h > 0)
{
	print "\n\n Script: 01_phyloTAGs_align.pl \n Author: Alejandro Caro-Quintero \n Version: 2.0 \n\n Description:\n";
	print " -in (File with sequences for primer design in fasta format)  \n";
	print " -no_distance (option for alignment by codons without the extraction by identity) \n";
	print " -distance (file with the pairwise realtedness distance)\n";
	print " -reference (identifier of reference sequence, e.g., NC_012345) \n";
	print " -idlow (lower value for pairwise distance to be evaluated) \n";	
	print " -idhigh (higher value for pairwise distance to be evaluated) \n\n";
}
else
{
    
    if($nodist < 1)
    {

	my%ident;
	my$low=$lowest;
	my$high=$low+1;
	my%hash;
	if ($low < 1)
	{
	    print "\n\n   **** Please specify all parameter (see parameters -help) **** \n\n";
	}	
	else	
	{
	    $pass=1;
	}
	## LOADING IDS FROM ALIGNMENT TO ENSURE THE EXISTENCE OF SEQ IN PROTEINS ALIGNMNET 
	%ident=loading($alignment);
  	$/="\n";
	if($pass > 0)
	{
	    open(FILE1, $SrRNA_id)||die("Could not open file $SrRNA_id 16S\n");
	    while(<FILE1>)
	    {
		chomp;
		my$line=$_;
		my($o1, $o2, $id)=split("\t", $line);
		if($o1 eq $reference || $o2 eq $reference)
		{			
		    if($id <= $highest && $id >= $lowest)
		    {
			if(exists($ident{$o1}) && exists($ident{$o2}))
			{
			    $hash{$o1}=$ident{$o1};
			    $hash{$o2}=$ident{$o2};
			}
			#else{ print "$o1\t$o2\n";}
		    }
		}
	    }
	    close(FILE1);	
	    
	    ##SAVING GENES RELATED TO REFERENCE
	    
	    open(OUTPUT,">>filtered\_$alignment");
	    {
		foreach my$guo (keys(%hash))
		{
		    if(exists($hash{$guo}))
		    {
			my$fsa=makefasta($guo,$hash{$guo});
			print OUTPUT "$fsa\n";
		    }
		}	
	    }
	    close(OUTPUT);
	}
		
	
	translate("filtered\_$alignment", "temp\_pro");
	system ("$clustalw -INFILE=temp\_pro -OUTPUT=FASTA -OUTORDER=INPUT -OUTFILE=temp\_pro\_aln");
	system ("$pal2nal temp\_pro\_aln filtered\_$alignment -output fasta >align_by_codons\_$alignment");
	system ("more align_by_codons\_$alignment \| tr \"_\" \" \"|awk \'\{if(\$1\~ \/>\/) print \$1\"_\"\$2;else print}\'>aligned_by_codons\.fas");
        system ("rm temp\_pro temp\_pro\_aln align_by_codons\_$alignment");
	makeolt("aligned_by_codons\.fas","aligned_by_codons\.olt"); 
	undef %ident;
        } 

  ##############################################
  ####### END OF DOING ALIGNMENT BY CODONS #####  
  ##############################################
  

    else
    { 
    translate("$alignment", "temp\_pro");
    system ("$clustalw -INFILE=temp\_pro -OUTPUT=FASTA -OUTORDER=INPUT -OUTFILE=temp\_pro\_aln");
    system ("$pal2nal temp\_pro\_aln $alignment -output fasta >align_by_codons\_$alignment");
    system ("more align_by_codons\_$alignment \| tr \"_\" \" \"|awk \'\{if(\$1\~ \/>\/) print \$1\"_\"\$2;else print}\'>aligned_by_codons\.fas");
    system ("rm temp\_pro temp\_pro\_aln align_by_codons\_$alignment");
    makeolt("aligned_by_codons\.fas","aligned_by_codons\.olt"); 
    #undef %ident;
    }
}
#################################################
#SUBROUTINES#####################################
#################################################		
sub loading
  {
    my $file1=$_[0];
    my $n=0; 
    my $data;
    my %load;
    {
      #gets fragment between signs
      $/=">";
      open (FILE,$file1); 
      while(<FILE>)
	{
	  my$temper=$_;	
	  my (@input)=split("\n",$temper);
	  my $l=@input;
	  my $t=$l-1;
	  my $name_seq =$input[0];
	  my $seq=join('',@input[1..$t]);
	  $seq=~s/>//g;
	  $name_seq=~s/>//g;
	  if($name_seq=~/>/)
	    {
	    $load{$name_seq}=$seq;
	    }
	  else
	    {
	     $load{$name_seq}=$seq;
	     #print "$load{$name_seq}\n";
	     }
	}
      return %load;
      close (FILE);
    }
  }

  sub makefasta 
    {
      my$fasta="";
      my$id=$_[0];
      my$raw_seq=$_[1];
      my$i=0;
      while ($i < length($raw_seq))
	{
	  my$seq_fasta = substr($raw_seq, $i, 60);
	  chomp $seq_fasta;
	  $fasta = "$fasta$seq_fasta\n";
	  $i= $i+60;
	}
      my$final_seq=">$id\n$fasta";
      return $final_seq;
    }
    
    sub makeolt
    {
$/="\n";    	
	my$file1=$_[0];
    	my$out=$_[1];
    	my$f="";
    	open(OUT,">>$out");
    	open (FILE,$file1)||die("Could not open file $file1\n");; 
      	while(<FILE>)
     	{
   			chomp;
    		my$data2=$_;
    		if($data2=~/^>(.+)/)
      		{
				$data2=~ s/>//;
				if($f=~/[A-Z]/)
				{
				print OUT "$f\n";
				}
				if($data2=~/[A-Z]/)
				{
				print OUT "$data2\t";
				}
				$f="";
      		}
    		else
      		{
			$f=$f.$data2;
      		}
  		}
	print OUT "$f\n";
    }

sub translate
{
    my $file1=$_[0];
    my $file2=$_[1];
    my $n=0; 
    my $data;
    my %load;
    my %transla;
    {
#####################
#LOAD SEQ INTO A HASH
#####################
	$/=">";
	open(OUT, ">>$file2");
	open(FILE,$file1)||die("Could not open $file1\n"); 
	while(<FILE>)
	{
	    my$temper=$_;	
	    my (@input)=split("\n",$temper);
	    my $l=@input;
	    my $t=$l-1;
	    my $name_seq =$input[0];
	    my $seq=join('',@input[1..$t]);
	    $seq=~s/>//g;
	    $name_seq=~s/>//g;
	    if($name_seq=~/>/)
	    {
		$load{$name_seq}=$seq;
	    }
	    else
	    {
		$load{$name_seq}=$seq;
	    }
	}	
	close (FILE);
    }
#########################
#TRANSLATE TO AMINO ACIDS
#########################
    foreach my$k (keys(%load))
    {
	my $heads=$k;    	
	my $sequ=$load{$k};
	my $tran="";
    	my%aacode = 
	    (
	     TTT => "F", TTC => "F", TTA => "L", TTG => "L",
	     TCT => "S", TCC => "S", TCA => "S", TCG => "S",
	     TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
	     TGT => "C", TGC => "C", TGA => "*", TGG => "W",
	     CTT => "L", CTC => "L", CTA => "L", CTG => "L",
	     CCT => "P", CCC => "P", CCA => "P", CCG => "P",
	     CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
	     CGT => "R", CGC => "R", CGA => "R", CGG => "R",
	     ATT => "I", ATC => "I", ATA => "I", ATG => "M",
	     ACT => "T", ACC => "T", ACA => "T", ACG => "T",
	     AAT => "N", AAC => "N", AAA => "K", AAG => "K",
	     AGT => "S", AGC => "S", AGA => "R", AGG => "R",
	     GTT => "V", GTC => "V", GTA => "V", GTG => "V",
	     GCT => "A", GCC => "A", GCA => "A", GCG => "A",
	     GAT => "D", GAC => "D", GAA => "E", GAG => "E",
	     GGT => "G", GGC => "G", GGA => "G", GGG => "G", 
	     "---" => "",
	    );
	
	$sequ=~ s/(...)/$1\_/g;
	my@codons=split("_",$sequ);
	foreach my$h (@codons)
	{
	    $tran="$tran"."$aacode{$h}";
	}
#################      	
#CONVERT TO FASTA
#################	
	my$fasta="";
      	my$i=0;
      	while ($i < length($tran))
	{
	    my$seq_fasta = substr($tran, $i, 60);
	    chomp $seq_fasta;
	    $fasta = "$fasta$seq_fasta\n";
	    $i= $i+60;
	}
	if($k=~ /[A-Z]|[0-9]/)
		{	
		my$final_seq=">$k\n$fasta";
		print OUT "$final_seq\n";
		}
    	}
 }
    
