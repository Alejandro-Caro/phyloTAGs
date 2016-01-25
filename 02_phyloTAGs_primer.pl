#!/usr/bin/perl -w
use strict;
use Getopt::Long;

use POSIX qw(strftime);
my $date = strftime "%m/%d/%Y", localtime;


################################################################################################################
# Script Name:02_phyloTAGs_primer.pl
# Script Version: 1.0
# Author: A. Caro-Quintero
# Date: July 2015
#
# Description: This script designs all posible primers forward and reverse with the required degeneracies from a 
# multi "fasta" or multi "olt" file, it also calculates the number of degeneracies.
#
# Parameters: 
# --in
# --format 
# --codons
# --slide
# --consensus
# 
##################################################################################################################

#Options variables

my$in="";#list with alignments per codon (output from "" olt)
my$co=3;
my$sl=0;
my$h=0;
my$j=0;
my$format="olt";
my$gredy=0;
my$consensus=100;
my$remove=0;

#other variables
my$control=0;
my$control2=0;
my$coor1=0;
my$coor2=0+($co*3);
my%hash;
my%temp;
my$long;
my$n=0;
my%place;
my%secuencia;
my%sec_dege;
my$range_low=0;
my$range_high=0;
my%order;
my%dege_code=("A", "A", "C", "C", "T","T","G","G", "U","U",  "ATCG","N",	"ATGC","N",	"ACTG","N",	"ACGT","N",	"AGTC","N",	
"AGCT","N",	"TAGC","N",	"TACG","N",	"TCGA","N",	"TCAG","N",	
"TGCA","N",	"TGAC","N",	"CATG","N",	"CAGT","N",	"CTAG","N",	
"CTGA","N",	"CGAT","N",	"CGTA","N",	"GACT","N",	"GATC","N",	
"GTCA","N",	"GTAC","N",	"GCTA","N",	"GCAT","N",	
"AG", "R", "GA", "R",
"TC", "Y", "CT", "Y",
"AC", "M", "CA", "M", 
"CG", "S", "GC", "S",
"AT", "W", "TA", "W",
"TG", "K", "GT", "K", 
"ACG","V",	"AGC","V",	"CAG","V",	"CGA","V",	"GAC","V",	"GCA","V",
"ATG","D",	"AGT","D",	"TAG","D",	"TGA","D",	"GAT","D",	"GTA","D",	 "ATC","H",	"ACT","H",	"TAC","H",	"TCA","H",	"CAT","H",	"CTA","H",
"TCG","B",	"TGC","B",	"CTG","B",	"CGT","B",	"GTC","B",	"GCT","B",
);
my%complement=("A", "T", "C", "G", "T","A","G","C", "U","A", "Y","R","R","Y","S","S","W","W","K","M","M","K","B","V","V","B","D","H","H","D","N","N", "-", "-");




# if not specified the values will be taken by default 
GetOptions("in:s" => \$in, "format:s" => \$format, "codons:s" => \$co, "slide:s" => \$sl,"consensus:s"=>\$consensus, "greedy"=>\$gredy, "help"=>\$h);

# printing help 
if ($h > 0)
{
	print "\n\n Script:  02_phyloTAGs_primer.pl \n Author: Alejandro Caro-Quintero \n Version: 1.0 \n\n Description:\n -in (File with mutiple sequence alignment) \n -format (file in \"fasta\" or \"olt\" format, olt by default) \n -codons (primer size in codons, e.g., 7 codons=7*3 =21 nucleotides)  \n -slide (sliding window size in bp) \n -consensus (percentage for calculation of consensus nucleotides) \n -greedy (explores all windows of codon size)\n\n";


}
else
{

	if($in=~ /[a-z]/ )
	{
		if (-e $in) {$control=1;}
		else{print "file $in not found\n";}
	}	
	else
	{
	print "Please specified input file (e.g --in filename)\n";
	}
}

open(OUT,">>phyloTAGs_primers_table.txt");
my$header="##-- phyloTAGs run on $date --  \n##  Authors: Caro-Quintero A. and Ochman H. 2015 \n\n"."Position_window_(bp)"."\t"."Degeneracies"."\t"."Bases_per_position"."\t"."Forward_primer"."\t"."Reverse_primer"."\n";
print OUT "$header";


if($format=~ /fasta/)
	{
	my$infile="$in";
	$in=~ s/fas/_olt/g;
	makeolt("$infile","$in");
	$control2=1;
	}
elsif($format=~ /olt/)
	{
	$control2++;
	}
	else
	{
	print "Unknow file format $in \n";	
	$control=0;
	}

if($gredy < 1 && $control2 > 0)

{
	#for nucleotide pr = 0
		
	if ($control > 0 ) 
	{
	$remove=1;
		open(FILE,$in)||die ("*****Could not find file $in *****\n");
		while(<FILE>)
		{
			chomp;
			my$l=$_;
			my($h,$s)=split("\t",$l);
			my@t=split("",$s);
			my$o=@t;
			$o=$o-1;
			my@slt=@t[$sl..$o];
			my$w=0;
			while ($sl>$w)
			{
				push(@slt,"-");
				$w++;
			}	
			my$kl=(join"", @slt);
			$temp{$h}=$kl;
		}
		open (OUT1, ">$sl\_$in");
		foreach my$x (keys(%temp))
		{
			print OUT1 "$x\t$temp{$x}\n";	
		}
		close(OUT1);
	}

	if ($control > 0 && $control2 > 0 ) 
	{
		open(FILE,"$sl\_$in")||die ("*****Could not find file $in *****\n");
		while(<FILE>)
		{
			chomp;
			my$line=$_;
			#print "$line\n";
			my($head,$seq)=split("\t",$line);
			$seq=~ s/(...)/$1\_/g;
			my@temp=split("_",$seq);
			$long=@temp;
			print "$long\n";
			$hash{$head}=\@temp;
		}	
		while($n<$long)
		{
			my%pos1;
			my%pos2;
			my%pos3;
			my$gap=1;
			my$seq="(";
			my$seq1="";
			my$seq2="";	
			my$seq3="";
			my$d_g="";
			my@numb=(keys(%hash));
			my$numb_seqs=(@numb);
			print "$numb_seqs\n";
			foreach my $keys (keys(%hash))
			{
				my$tempo=@{$hash{$keys}}[$n];
				#print "$tempo\n";
				my@temp=split("",$tempo);
				my$f=@temp;
				if ($f >2)
				{
					#counting number of times base is find per codon position 1
					if(exists($pos1{$temp[0]}))
					{
						$pos1{$temp[0]}=$pos1{$temp[0]}+1;
					}
					else 
					{
					$pos1{$temp[0]}=1;
					}
					#counting number of times base is find per codon position 2					
					
					if(exists($pos2{$temp[1]}))
					{
						$pos2{$temp[1]}=$pos2{$temp[1]}+1;
					}
					else 
					{
					$pos2{$temp[1]}=1;
					}
					#counting number of times base is find per codon position 3					
					if(exists($pos3{$temp[2]}))
					{
						$pos3{$temp[2]}=$pos3{$temp[2]}+1;
					}
					else 
					{
					$pos3{$temp[2]}=1;
					}
					
				}
			}
			#calculating consensus using defined percentage 
			foreach my$pa (keys(%pos1))
				{
				 if (($pos1{$pa}/$numb_seqs)*100 >$consensus)
					{
					 $pos1{"consensus"}=$pa;
					}
						
				}

			foreach my$pb(keys(%pos2))
				{
					if (($pos2{$pb}/$numb_seqs)*100 > $consensus)
					{
					 $pos2{"consensus"}=$pb;
					}
						
				}

			foreach my$pc (keys(%pos3))
				{
				 if (($pos3{$pc}/$numb_seqs)*100 > $consensus)
					{
					 $pos3{"consensus"}=$pc;
					}
						
				}
			#IF CONSENSUS IS IN HASH DELETE THE OTHER KEYS AND SUBSTITUTE BY VALUE OF CONSENSUS 
			#replacing for consensus when necessary	
			if (exists($pos1{"consensus"}))
			{
				my$temp= $pos1{"consensus"};
				%pos1=();
				$pos1{$temp}="";
			}
			
			if (exists($pos2{"consensus"}))
			{
				my$temp= $pos2{"consensus"};
				%pos2=();
				$pos2{$temp}="";
			}
			
			if (exists($pos3{"consensus"}))
			{
				my$temp= $pos3{"consensus"};
				%pos3=();
				$pos3{$temp}="";
			}
			
			my$p1=(keys %pos1);
			my$p2=(keys %pos2);	
			my$p3=(keys %pos3);
	
		
			### POS1 ###
			foreach my$h (keys %pos1)
			{
				if ($h eq "-"){$gap=0;}
				$seq="$seq"."$h";
				$seq1="$seq1"."$h";
			}
			##
			if(exists($dege_code{$seq1}))
			{
				$d_g=$d_g."$dege_code{$seq1}";
			}	
			else {$d_g=$d_g."-";}
			$seq=$seq.")(";
		
			### POS2 ###
		
			foreach my$u (keys %pos2)
			{
				if ($u eq "-"){$gap=0;}
				$seq="$seq"."$u";
				$seq2="$seq2"."$u";
			}
			$seq=$seq.")(";	
		
			##
			if(exists($dege_code{$seq2}))
			{
				$d_g=$d_g."$dege_code{$seq2}";
			}
			else {$d_g=$d_g."-";}
		
		
			### POS3 ###
		
		
			foreach my$o (keys %pos3)
			{
			if ($o eq "-"){$gap=0;}
			$seq="$seq"."$o";
			$seq3="$seq3"."$o";
			}
			$seq=$seq.")";

		
			### 
			if(exists($dege_code{$seq3}))
			{
				$d_g=$d_g."$dege_code{$seq3}";
			}
			else {$d_g=$d_g."-";}
		
		
			$secuencia{$n}="$seq";
			$sec_dege{$n}="$d_g";
		
			
			my$c=$p1*$p2*$p3*$gap;
				
			$place{$n}=$c;
			$n++;
			undef(%pos1);
			undef(%pos2);
			undef(%pos3);
		}
	


		$range_high=$range_low+$co-1;

		my@y=sort{$a<=>$b}(keys(%place));
		while($range_high< $long)	
		{
	
			my@sub=@y[$range_low..$range_high];
			my$tot=1;
			my$rest="";
			my$peace="";
			my@piece;
			foreach  my$key (@sub)
			{
				$tot=$tot*$place{$key};
				$rest=$rest."$secuencia{$key}";
				$peace=$peace."$sec_dege{$key}";
		
			}
			my$r1=3*($range_low+1)-2;
			my$r2=3*($range_high+1);
			#print "$peace\n";
			my@splt=split("",$peace);
			foreach my$kl (@splt)
			{
				#print "$kl\n";
				push(@piece,"$complement{$kl}");
			}
			my@reverse_piece=reverse(@piece);
			my$rev_comp=join("",@reverse_piece);
		
			my$h1=$r1+$sl;
			my$h2=$r2+$sl;
			$order{$h1}="$h1"."_"."$h2\t$tot\t$rest\t$peace\t$rev_comp\n";
			$range_high=$range_high+$co;
			$range_low=$range_low+$co;
		}		

	}
if($remove>0)	{system("rm \*_aligned_by_codons\*");
		}
}
#if gredy is on DO
else 
{
		
	while ($sl<($co*3))
	{
		#for nucleotide pr = 0
	#	print "$sl\n";
		if ($control > 0 && $control2 > 0 ) 
		{
			$remove=1;
			open(FILE,$in)||die ("*****Could not find file $in *****/n");
			while(<FILE>)
			{
				chomp;
				my$l=$_;
				my($h,$s)=split("\t",$l);
				my@t=split("",$s);
				my$o=@t;
				$o=$o-1;
				my@slt=@t[$sl..$o];
				my$w=0;
				while ($sl>$w)
				{
					push(@slt,"-");
					$w++;
				}	
				my$kl=(join"", @slt);
				$temp{$h}=$kl;
			}
			open (OUT1, ">$sl\_$in");
			foreach my$x (keys(%temp))
			{
				print OUT1 "$x\t$temp{$x}\n";	
			}
			close(OUT1);


		}
	$n=0;
	$range_low=0;
	$range_high=0;

		#for nucleotide pr = 0
		if ($control > 0 && $control2 > 0 ) 
		{
			open(FILE,"$sl\_$in")||die ("*****Could not find file $in *****/n");
			while(<FILE>)
			{	
				chomp;
				my$line=$_;
				#print "$line\n";
				my($head,$seq)=split("\t",$line);
				$seq=~ s/(...)/$1\_/g;
				my@temp=split("_",$seq);
				$long=@temp;
			#	print "$long\n";
				$hash{$head}=\@temp;
			}	
	
			while($n<$long)
			{
				my%pos1;
				my%pos2;
				my%pos3;
				my$gap=1;
				my$seq="(";
		
				my$seq1="";
				my$seq2="";	
				my$seq3="";
				my$d_g="";
				my@numb=(keys(%hash));
				my$numb_seqs=(@numb);
				foreach my $keys (keys(%hash))
				{
					my$tempo=@{$hash{$keys}}[$n];
					#print "$tempo\n";
					my@temp=split("",$tempo);
					my$f=@temp;
					if ($f >2)
					{
						#counting number of times base is find per codon position 1
						if(exists($pos1{$temp[0]}))
						{
							$pos1{$temp[0]}=$pos1{$temp[0]}+1;
						}	
						else 
						{
							$pos1{$temp[0]}=1;
						}
						#counting number of times base is find per codon position 2					
					
						if(exists($pos2{$temp[1]}))
						{
							$pos2{$temp[1]}=$pos2{$temp[1]}+1;
						}
						else 
						{
							$pos2{$temp[1]}=1;
						}
						#counting number of times base is find per codon position 3					
						if(exists($pos3{$temp[2]}))
						{
							$pos3{$temp[2]}=$pos3{$temp[2]}+1;
						}
						else 
						{
							$pos3{$temp[2]}=1;
						}
					}
				}
				
					#calculating consensus using defined percentage 
			foreach my$pa (keys(%pos1))
				{
				 if (($pos1{$pa}/$numb_seqs)*100 >$consensus)
					{
					 $pos1{"consensus"}=$pa;
					}
						
				}

			foreach my$pb(keys(%pos2))
				{
					if (($pos2{$pb}/$numb_seqs)*100 > $consensus)
					{
					 $pos2{"consensus"}=$pb;
					}
						
				}

			foreach my$pc (keys(%pos3))
				{
				 if (($pos3{$pc}/$numb_seqs)*100 > $consensus)
					{
					 $pos3{"consensus"}=$pc;
					}
						
				}
			#IF CONSENSUS IS IN HASH DELETE THE OTHER KEYS AND SUBSTITUTE BY VALUE OF CONSENSUS 
			#replacing for consensus when necessary	
			if (exists($pos1{"consensus"}))
			{
				my$temp= $pos1{"consensus"};
				%pos1=();
				$pos1{$temp}="";
			}
			
			if (exists($pos2{"consensus"}))
			{
				my$temp= $pos2{"consensus"};
				%pos2=();
				$pos2{$temp}="";
			}
			
			if (exists($pos3{"consensus"}))
			{
				my$temp= $pos3{"consensus"};
				%pos3=();
				$pos3{$temp}="";
			}
									
				my$p1=(keys %pos1);
				my$p2=(keys %pos2);	
				my$p3=(keys %pos3);
	
		
				### POS1 ###
				foreach my$h (keys %pos1)
				{
					if ($h eq "-"){$gap=0;}
					$seq="$seq"."$h";
					$seq1="$seq1"."$h";
				}
				##
				if(exists($dege_code{$seq1}))
				{
					$d_g=$d_g."$dege_code{$seq1}";
				}
				else {$d_g=$d_g."-";}
		
				$seq=$seq.")(";
		
				### POS2 ###
		
				foreach my$u (keys %pos2)
				{
					if ($u eq "-"){$gap=0;}
					$seq="$seq"."$u";
					$seq2="$seq2"."$u";
				}
				$seq=$seq.")(";	
		
				##
				if(exists($dege_code{$seq2}))
				{
					$d_g=$d_g."$dege_code{$seq2}";
				}
				else {$d_g=$d_g."-";}
		
		
				### POS3 ###
		
		
				foreach my$o (keys %pos3)
				{
					if ($o eq "-"){$gap=0;}
					$seq="$seq"."$o";
					$seq3="$seq3"."$o";
				}
				$seq=$seq.")";
				### 
				if(exists($dege_code{$seq3}))
				{
					$d_g=$d_g."$dege_code{$seq3}";
				}
				else {$d_g=$d_g."-";}
		
				$secuencia{$n}="$seq";
				$sec_dege{$n}="$d_g";
		
			
				my$c=$p1*$p2*$p3*$gap;
				
				$place{$n}=$c;
				$n++;
				undef(%pos1);
				undef(%pos2);
				undef(%pos3);
			}
	


			$range_high=$range_low+$co-1;

			my@y=sort{$a<=>$b}(keys(%place));
			while($range_high< $long)	
			{
	
				my@sub=@y[$range_low..$range_high];
				my$tot=1;
				my$rest="";
				my$peace="";
				my@piece;
				foreach  my$key (@sub)
				{
					$tot=$tot*$place{$key};
					$rest=$rest."$secuencia{$key}";
					$peace=$peace."$sec_dege{$key}";
		
				}
				my$r1=3*($range_low+1)-2;
				my$r2=3*($range_high+1);
				#print "$peace\n";
				my@splt=split("",$peace);
				foreach my$kl (@splt)
				{
					#print "$kl\n";
					push(@piece,"$complement{$kl}");
				}
				my@reverse_piece=reverse(@piece);
				my$rev_comp=join("",@reverse_piece);
		
				my$h1=$r1+$sl;
				my$h2=$r2+$sl;
				$order{$h1}="$h1"."_"."$h2\t$tot\t$rest\t$peace\t$rev_comp\n";
				$range_high=$range_high+$co;
				$range_low=$range_low+$co;
			}		

		}
	$sl++;
	}
if($remove>0)	{
		system("rm \*_aligned_by_codons\*");
		}
}	
foreach my$keyo (sort { $a <=> $b }(keys(%order)))
	{
	my$ho=$order{$keyo};
	print OUT "$ho";
	}
	#system("rm \*aligned_by_codons\*");


sub makeolt
    {
   	
	my$file1=$_[0];
    	my$out=$_[1];
    	my$f="";
	my$cc=0;
	my$data2="";
    	open(OUT2,">>$out");
    	open (FILE,$file1)||die("Could not open file $file1\n");; 
	while(<FILE>)      	  
	{
  		chomp;
    		$data2=$_;
    		if($data2=~/^>(.+)/)
      		{
			$data2=~ s/>//;
			if($cc>0)
			{
			print OUT2 "$f\n";
			print OUT2 "$data2\t";
			$f="";
			}
			else
			{
			print OUT2 "$f";
			print OUT2 "$data2\t";
			$f="";
			$cc++;
			}
			
			
      		}
    		else
      		{
			$f=$f.$data2;
      		}		
  	}
	print OUT2 "$f\n";
	#return $cc;
	close(OUT2);
	close (FILE);
    }
