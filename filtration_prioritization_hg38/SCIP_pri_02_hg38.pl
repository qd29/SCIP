#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Cwd 'abs_path';
my %opts;
getopt ('n:',\%opts);
my $name=$opts{"n"};

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

open out1, ">$dir/user_data/$name.pipeline_summary.txt";
open out2, ">$dir/user_data/$name.pipeline_rerun.txt";
my %dedup;
open file1, "<$dir/user_data/$name.variant_list.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 my $temp="$split1[0].$split1[1].$split1[2].$split1[3].$split1[4]";
 unless (exists $dedup{$temp}){
  $dedup{$temp}=1;

  if (-e "$dir/d1stat/$temp.txt"){
   open file2, "<$dir/d1stat/$temp.txt";
   my $line=<file2>; chomp $line; close file2;
   my @split2=split /\t/,$line;
   my $length=$split1[3]-$split1[2]+1;

   my $class="Manual";
   if ($split2[2] ne "NA"){
    if ($split2[2]<=0.1 && $split1[4] eq "DEL"){
     $class="Manual-HomDEL (DEL,depth_ratio<=0.1)";
    }
    elsif ($length>=2000 && $split2[3]<=40 && $split2[6]==0 && $split2[7]==0 && $split2[9]==0){
     $class="Failed-1 (length>=2kb,MQ<=40,no_supporting_or_opposite_PE,no_SE)";
    }
    elsif ($split2[3]>=50 && $split2[6]>=2 && $split1[4] eq "DEL" && $split2[2]<=0.75){
     $class="Passed-1 (DEL,depth_ratio<=0.75,MQ>=50,supporting_PE>=2)";
    }
    elsif ($split2[3]>=50 && $split2[6]>=2 && $split1[4] eq "DUP" && $split2[2]>=1.25){
     $class="Passed-2 (DUP,depth_ratio>=1.25,MQ>=50,supporting_PE>=2)";
    }
    elsif ($split2[3]>=50 && $split2[9]>=2 && $split1[4] eq "DEL" && $split2[2]<=0.75){
     $class="Passed-3 (DEL,depth_ratio<=0.75,MQ>=50,supporting_SE>=2)";
    }
    elsif ($split2[3]>=50 && $split2[9]>=2 && $split1[4] eq "DUP" && $split2[2]>=1.25){
     $class="Passed-4 (DUP,depth_ratio>=1.25,MQ>=50,supporting_SE>=2)";
    }
    elsif ($length>=10000 && $split2[3]>=50 && $split1[4] eq "DEL" && $split2[2]<=0.7){
     $class="Passed-5 (DEL,length>=10kb,depth_ratio<=0.7,MQ>=50)";
    }
    elsif ($length>=10000 && $split2[3]>=50 && $split1[4] eq "DUP" && $split2[2]>=1.3){
     $class="Passed-6 (DUP,length>=10kb,depth_ratio>=1.3,MQ>=50)";
    }
   }
   
   if ($split2[7]>=2){
    $class="COMPLEX-$class";
   }

   open file2, "gunzip -c $dir/app_temp_file/$split1[0]/$temp.script08_file1.txt.gz |";
   my $line2=<file2>; chomp $line2; close file2;
   my @split3=split /\t/,$line2;
   unless ($split3[4]){
    $split3[4]="";
   }
   
   unless ($split3[4]=~/FULL_gnomAD_SV_common_var_OVERLAP/){
    my $priority=$split3[1]+100*$split3[2];
    print out1 "$class\t$_\t$line\t$temp\t$priority\t$split3[2]\t$split3[3]\n";
   }
  }
  else{
   print out2 "$_\n";
  }
 }
}
close file1;
close out1;
close out2;

exit 2;
