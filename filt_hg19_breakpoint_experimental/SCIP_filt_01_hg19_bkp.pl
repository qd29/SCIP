#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my %opts;
getopt ('n:',\%opts);
my $name=$opts{"n"};

my @type=("DEL","DUP","INS","BND");
my @chr;
for (my $i=1; $i<=22; $i++){
 push @chr, $i;
}
push @chr, "X";
# parallelization possible if this script is run by chromosome and/or variant type

open out1, ">./user_data/$name.bkp.filt_step01.txt";
foreach my $chr (@chr){
 foreach my $type (@type){
  my @reg;
  open file1, "<./hg19_files/hg19_recurrent_breakpoint.txt";
  while (<file1>){
   chomp;
   my @split1=split /\t/,$_;
   if ($split1[6] eq $type && $split1[0] eq $chr){
    my $temp="$split1[1]|$split1[2]|$split1[3]|$split1[4]";
    push @reg, $temp;
   }
  }
  close file1;

  open file1, "<./user_data/$name.bkp.unfiltered_CNV.txt";
  while (<file1>){
   chomp;
   my @split1=split /\t/,$_;
   my $st=0;
   if ($split1[0] eq $chr && $split1[3] eq $type){
    LOOP1: foreach my $reg (@reg){
     my @split2=split /\|/,$reg;
     my $exp1=$split2[0]-5000; my $exp2=$split2[1]+5000;
     my $exp3=$split2[2]-5000; my $exp4=$split2[3]+5000;
     if ($split1[1]>=$exp1 && $split1[1]<=$exp2 && $split1[2]>=$exp3 && $split1[2]<=$exp4){
      $st=1;
     }
     if ($split1[1]>=$split2[0] && $split1[1]<=$split2[1] && $split1[2]>=$split2[2] && $split1[2]<=$split2[3]){
      $st=2; last LOOP1;
     }
    }

    if ($st<2){
     print out1 "$_\t$st\n";
    }
   }
  }
  close file1;
 }
 print "SCIP Filtration Module script 01 processing hg19 (breakpoint-based) chr$chr\n";
}
close out1;

exit 2;
