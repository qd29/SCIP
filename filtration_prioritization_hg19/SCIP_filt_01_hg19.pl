#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my %opts;
getopt ('n:',\%opts);
my $name=$opts{"n"};
my $ext=0; # extend percentage

my @chr;
for (my $i=1; $i<=22; $i++){
 push @chr, $i;
}
push @chr, "X";
# parallelization possible if this script is run by chromosome

open out1, ">./user_data/$name.filt_step01.txt";
foreach my $chr (@chr){
my @shared;
open file1, "gunzip -c ./hg19_files/hg19_gap.txt.gz |";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[1] eq "chr$chr"){
  my $temp="$split1[2]|$split1[3]";
  push @shared, $temp;
 }
}
close file1;

open file1, "gunzip -c ./hg19_files/hg19_genomicSuperDups.txt.gz |";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[1] eq "chr$chr"){
  my $temp="$split1[2]|$split1[3]";
  push @shared, $temp;
 }
 if ($split1[7] eq "chr$chr"){
  my $temp="$split1[8]|$split1[9]";
  push @shared, $temp;
 }
}
close file1;

print "SCIP Filtration Module script 01 processing hg19 chr$chr\n";
my (@del,@dup);
open file1, "<./hg19_files/hg19_recurrence_region.DEL.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq $chr){
  my $len=$split1[2]-$split1[1]+1;
  $split1[1]-=int($ext*$len); $split1[2]+=int($ext*$len);
  push @del, "$split1[1]|$split1[2]";
 }
}
close file1;

open file1, "<./hg19_files/hg19_recurrence_region.DUP.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq $chr){
  my $len=$split1[2]-$split1[1]+1;
  $split1[1]-=int($ext*$len); $split1[2]+=int($ext*$len);
  push @dup, "$split1[1]|$split1[2]";
 }
}
close file1;

open file1, "<./user_data/$name.unfiltered_CNV.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq $chr){
  my $out=1;
  if ($split1[3] eq "DEL"){
   LOOP2: foreach my $shared (@shared){
    my @split2=split /\|/,$shared;
    if ($split2[0]<=$split1[1] && $split2[1]>=$split1[2]){
     $out=0; last LOOP2;
    }
   }
 
   if ($out==1){
    LOOP1: foreach my $del (@del){
     my @split2=split /\|/,$del;
     if ($split2[0]<=$split1[1] && $split2[1]>=$split1[2]){
      $out=0; last LOOP1;
     }
    }
   }
  }
 
  if ($split1[3] eq "DUP"){
   LOOP2: foreach my $shared (@shared){
    my @split2=split /\|/,$shared;
    if ($split2[0]<=$split1[1] && $split2[1]>=$split1[2]){
     $out=0; last LOOP2;
    }
   }
 
   if ($out==1){
    LOOP1: foreach my $dup (@dup){
     my @split2=split /\|/,$dup;
     if ($split2[0]<=$split1[1] && $split2[1]>=$split1[2]){
      $out=0; last LOOP1;
     }
    }
   }
  }
 
  if ($out==1){
   print out1 "$_\n";
  }
 }
}
close file1;
}
close out1;

system ("perl ./SCIP_filt_02_hg19.pl -n $name");
system ("perl ./SCIP_filt_03_hg19.pl -n $name");
system ("perl ./SCIP_filt_04_hg19.pl -n $name");
my $temp_name="$name.hg19";
system ("perl ./SCIP_pri_01_hg19.pl -n $temp_name");
system ("perl ./SCIP_pri_02_hg19.pl -n $temp_name");

exit 2;
