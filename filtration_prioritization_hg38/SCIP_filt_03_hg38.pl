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

open out1, ">./user_data/$name.filtered.txt";
foreach my $chr (@chr){
 print "SCIP Filtration Module script 03 processing hg38 chr$chr\n";
 my (@del,@dup);
 open file1, "<./hg38_files/gnomADv2.1_commonCNV_hg38.gene.txt";
 while (<file1>){
  chomp;
  my @split1=split /\t/,$_;
  if ($split1[0] eq $chr){
   my $temp="$split1[1]|$split1[2]";
   if ($split1[3] eq "DEL"){
    push @del, $temp;
   }
   elsif ($split1[3] eq "DUP"){
    push @dup, $temp;
   }
  }
 }
 close file1;

 open file1, "<./user_data/$name.filt_step02.txt";
 while (<file1>){
  chomp;
  my @split1=split /\t/,$_;
  if ($split1[0] eq $chr){
   if ($split1[3] eq "DEL"){
    my $st=0;
 
    LOOP1: foreach my $del (@del){
     my @split2=split /\|/,$del;
     if ($split2[0]<=$split1[1] && $split2[1]>=$split1[2]){
      $st=1; last LOOP1;
     }
    }
 
    if ($st==0){
     print out1 "$_\n";
    }
   }
   elsif ($split1[3] eq "DUP"){
    my $st=0;
  
    LOOP1: foreach my $dup (@dup){
     my @split2=split /\|/,$dup;
     if ($split2[0]<=$split1[1] && $split2[1]>=$split1[2]){
      $st=1; last LOOP1;
     }
    }

    if ($st==0){
     print out1 "$_\n";
    }
   }
  }
 }
 close file1;
}
close out1;

exit 2;
