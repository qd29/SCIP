#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Cwd 'abs_path';
my %opts;
getopt ('n:u:',\%opts);
my $name=$opts{"n"};
my $num=$opts{"u"};
my $start_num=($num-1)*50000+1;
my $end_num=$num*50000;
# adjust the above number (50000) lower to achieve parallelization

my $cwd=abs_path();
my $dir;
open file1, "<$cwd/pipeline_config_hg38.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq "CURRENT_PATH"){
  $dir=$split1[1];
 }
}
close file1;

open out1, ">$dir/$name.$num.variant_list.txt";
my $ct=1;
open file1, "<$dir/user_data/$name.variant_list.txt";
while (<file1>){
 chomp;
 if ($ct>=$start_num && $ct<=$end_num){
  print out1 "$_\n";
 }
 $ct++;
}
close file1;
close out1;

open file1, "<$dir/$name.$num.variant_list.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 my $start=localtime();
 my $start_time=time();
 print "$start\t$_\n";

 # comment out this bracket if want to re-generate temporary files for all existing variants
 unless (-e "$dir/app_temp_file/$split1[0]/$split1[0].$split1[1].$split1[2].$split1[3].$split1[4].log.txt.gz"){
  system ("mkdir -p $dir/app_temp_file/$split1[0]"); 
  system ("perl $dir/SCIP_pri_03_hg38.pl -c $split1[1] -s $split1[2] -e $split1[3] -p $split1[0] -t $split1[4]");
  system ("perl $dir/SCIP_pri_04_hg38.pl -c $split1[1] -s $split1[2] -e $split1[3] -p $split1[0] -t $split1[4]");
  system ("perl $dir/SCIP_pri_05_hg38.pl -c $split1[1] -s $split1[2] -e $split1[3] -p $split1[0] -t $split1[4]");
  system ("perl $dir/SCIP_pri_06_hg38.pl -c $split1[1] -s $split1[2] -e $split1[3] -p $split1[0] -t $split1[4]");
  system ("perl $dir/SCIP_pri_07_hg38.pl -c $split1[1] -s $split1[2] -e $split1[3] -p $split1[0] -t $split1[4]");
  system ("perl $dir/SCIP_pri_08_hg38.pl -c $split1[1] -s $split1[2] -e $split1[3] -p $split1[0] -t $split1[4]");
  my $end=localtime();
  my $end_time=time();
  my $timediff=$end_time-$start_time;

  open out1, ">$dir/app_temp_file/$split1[0]/$split1[0].$split1[1].$split1[2].$split1[3].$split1[4].log.txt";
  print out1 "Start\t$start\nComplete\t$end\nDuration (seconds)\t$timediff\n";
  open file2, "<$dir/pipeline_config_hg38.txt";
  while (<file2>){
   chomp;
   my @split2=split /\t/,$_;
   unless ($split2[0] eq "REF_BAM" || $split2[0] eq "CURRENT_PATH" || $split2[0] eq "ALIGNMENT_PATH" || $split2[0] eq "PREFIX"){
    my @split4=split /\//,$split2[1];
    system ("md5sum $split2[1] > $dir/$split1[0].$split1[1].$split1[3].$split1[4].d2temp1");
    open file3, "<$dir/$split1[0].$split1[1].$split1[3].$split1[4].d2temp1";
    my $line=<file3>; chomp $line; my @split3=split /\s/,$line;
    close file3;
    system ("rm $dir/$split1[0].$split1[1].$split1[3].$split1[4].d2temp1");
    print out1 "$split2[0]\t$split4[$#split4]\t$split3[0]\n";
   }
  }
  close file2;
  close out1;
  system ("gzip -f $dir/app_temp_file/$split1[0]/$split1[0].$split1[1].$split1[2].$split1[3].$split1[4].log.txt");
 }
}
close file1;

system ("rm $dir/$name.$num.variant_list.txt");
exit 2;
