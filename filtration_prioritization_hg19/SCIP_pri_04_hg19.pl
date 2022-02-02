#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use Getopt::Std;
my %opts;
getopt ('p:c:s:e:t:',\%opts);
my $chr=$opts{'c'};
my $start=$opts{'s'};
my $end=$opts{'e'};
my $proband=$opts{'p'};
my $type=$opts{"t"};

my $start_original=$start;
my $end_original=$end;
my $exp_length=$end_original-$start_original+1;
my $ext=int(($end-$start+1)*0.5);
if ($ext<=1e5){
 $ext=1e5;
}
$start-=$ext; $end+=$ext;
if ($start<1){$start=1}

my $cwd=abs_path();
my ($dir,$refbam);
open file1, "<$cwd/pipeline_config.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq "CURRENT_PATH"){
  $dir=$split1[1];
 }
 if ($split1[0] eq "REF_BAM"){
  $refbam=$split1[1];
 }
}
close file1;

#unless (-e "$dir/app_temp_file/$proband/$proband.$chr.$start_original.$end_original.$type.script04_file1.txt.gz"){
 system ("samtools depth $refbam -r chr$chr:$start-$end > $dir/$proband.$chr.$start_original.$end_original.$type.d101temp1");
 system ("samtools view -F 0x400 $refbam chr$chr:$start-$end |cut -f1,2,3,4,5,6,7,8,9,12 -d\$'\t' > $dir/$proband.$chr.$start_original.$end_original.$type.d101temp2");

 my $min=1e20; my $max=0;
 open file1, "<$dir/$proband.$chr.$start_original.$end_original.$type.d101temp2";
 my (%mq,%mqtot);
 while (<file1>){
  chomp;
  my @split1=split /\t/,$_;
  my $read_end=$split1[3];

  if ($split1[5] ne "*"){
   my @split11=split /M/,$split1[5];
   my $ext=-1;
   foreach my $split11 (@split11){
    unless ($split11=~/[A-Za-z]$/){
     my @split12=split /[A-Za-z]/,$split11;
     $ext+=$split12[$#split12];
    }
   }
   $read_end=$split1[3]+$ext;
  }
 
  if ($split1[3]<=$min){$min=$split1[3]}
  if ($split1[3]>=$max){$max=$split1[3]}

  for (my $i=$split1[3]; $i<=$read_end; $i++){
   $mq{$i}+=$split1[4]; $mqtot{$i}++;
  }
 }
 close file1;
 
 open file1, "<$dir/$proband.$chr.$start_original.$end_original.$type.d101temp1";
 my %depth;
 while (<file1>){
  chomp;
  my @split1=split /\t/,$_;
 
  if ($split1[1]<=$min){$min=$split1[1]}
  if ($split1[1]>=$max){$max=$split1[1]}
 
  $depth{$split1[1]}=$split1[2];
 }
 close file1;
 
 open out1, ">$dir/app_temp_file/$proband/$proband.$chr.$start_original.$end_original.$type.script04_file1.txt";
 for (my $i=$min; $i<=$max; $i++){
  my $depth="NA"; my $mq="NA";
  if (exists $depth{$i}){
   $depth=$depth{$i};
  }
  if (exists $mq{$i}){
   $mq=int($mq{$i}*1e2/$mqtot{$i})/1e2;
  }
  print out1 "$i\t$depth\t$mq\n";
 }
 close out1;
 system ("gzip -f $dir/app_temp_file/$proband/$proband.$chr.$start_original.$end_original.$type.script04_file1.txt");
 system ("rm $dir/$proband.$chr.$start_original.$end_original.$type.d101temp1 $dir/$proband.$chr.$start_original.$end_original.$type.d101temp2");
#}

exit 2;
