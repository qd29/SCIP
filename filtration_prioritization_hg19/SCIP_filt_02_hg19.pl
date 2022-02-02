#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my %opts;
getopt ('n:',\%opts);
my $name=$opts{"n"};

my @chr;
for (my $i=1; $i<=22; $i++){
 push @chr, $i;
}
push @chr, "X";
# parallelization possible if this script is run by chromosome

open out1, ">./user_data/$name.filt_step02.txt";
foreach my $chr (@chr){
 print "SCIP Filtration Module script 02 processing hg19 chr$chr\n";
 my @exon;
 open file1, "<./hg19_files/hg19_coding_exons_v3_name_range.txt";
 while (<file1>){
  chomp;
  my @split1=split /\t/,$_;
  if ($split1[2] eq $chr){
   my @split2=split /\s/,$split1[5];
   foreach my $split2 (@split2){
    push @exon, $split2;
   }
  }
 }
 close file1;

 open file1, "<./user_data/$name.filt_step01.txt";
 while (<file1>){
  chomp;
  my @split1=split /\t/,$_;
  if ($split1[0] eq $chr){
   my $st=0;

   LOOP1: foreach my $exon (@exon){
    my @split2=split /\|/,$exon;
    if ($split2[2]>=$split1[1] && $split2[1]<=$split1[2]){
     $st=1; last LOOP1;
    }
   }

   if ($st==1){
    print out1 "$_\n";
   }
  }
 }
 close file1;
}
close out1;

exit 2;
