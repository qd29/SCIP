#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);
use Getopt::Std;
use Cwd 'abs_path';
my %opts;
getopt ('c:s:e:p:t:',\%opts);
my $chr=$opts{'c'};
my $start=$opts{'s'};
my $end=$opts{'e'};
my $proband=$opts{'p'};
my $type=$opts{"t"};

my $cwd=abs_path();
open file1, "<$cwd/pipeline_config.txt";
my ($gnomad_pext_file,$coding_exons_plotting,$dir);
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq "gnomAD_pext"){
  $gnomad_pext_file=$split1[1];
 }
 if ($split1[0] eq "coding_exons_plotting"){
  $coding_exons_plotting=$split1[1];
 }
 if ($split1[0] eq "CURRENT_PATH"){
  $dir=$split1[1];
 }
}

my @cyc=("d1temp_server");
my $preload=0;
my ($depth,$sam);
foreach my $cyc (@cyc){
 if (-e "$dir/$cyc/$proband.$chr.$start.$end.$type.SAM.txt.gz"){
  $preload=1;
  $depth="$dir/$cyc/$proband.$chr.$start.$end.$type.DEPTH.txt.gz";
  $sam="$dir/$cyc/$proband.$chr.$start.$end.$type.SAM.txt.gz";
 }
}
if ($preload==0){
 die "SAM file does not exist for $proband.$chr.$start.$end.$type\n";
}

my $min=1e20; my $max=0;
open file1, "gunzip -c $sam |";
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

open file1, "gunzip -c $depth |";
my %depth;
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;

 if ($split1[1]<=$min){$min=$split1[1]}
 if ($split1[1]>=$max){$max=$split1[1]}

 $depth{$split1[1]}=$split1[2];
}
close file1;

system ("tabix $gnomad_pext_file $chr:$min-$max > $dir/$proband.$chr.$start.$end.$type.d4temp1.txt");
open file1, "<$dir/$proband.$chr.$start.$end.$type.d4temp1.txt";
my %pext;
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 $pext{$split1[2]}=int($split1[3]*1e2)/1e2;
}
close file1;
system ("rm $dir/$proband.$chr.$start.$end.$type.d4temp1.txt");

open out1, ">$dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script05_file1.txt";
for (my $i=$min; $i<=$max; $i++){
 my $depth="NA"; my $mq="NA"; my $pext="NA";
 if (exists $depth{$i}){
  $depth=$depth{$i};
 }
 if (exists $pext{$i}){
  $pext=$pext{$i};
 }
 if (exists $mq{$i}){
  $mq=int($mq{$i}*1e2/$mqtot{$i})/1e2;
 }
 print out1 "$i\t$depth\t$mq\t$pext\n";
}
close out1;

my (%gene_range,%gene_level);
my $level=1;
open file1, "<$coding_exons_plotting";
open out1, ">$dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script05_file3.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[1] eq $chr && $split1[3]>=$min && $split1[2]<=$max){
  $gene_range{$split1[0]}="$split1[2]|$split1[3]";

  my @arr1;
  my @split4=split /\s/,$split1[4];
  foreach my $split4 (@split4){
   my @split5=split /\|/,$split4;
   push @arr1, $split5[1]; push @arr1, $split5[2];
  }
  @arr1=sort{$a<=>$b}@arr1;

  my $st1=0;
  my @gene_level=keys (%gene_level);
  @gene_level=shuffle(@gene_level);

  LOOP1: foreach my $gene_level (@gene_level){
   my $overlap=0;
   my @split2=split /\|/,$gene_level{$gene_level};
   foreach my $split2 (@split2){
    my @split3=split /\|/,$gene_range{$split2};
    if ($split1[3]>=$split3[0] && $split1[2]<=$split3[1]){
     $overlap=1;
    }
   }
   if ($overlap==0){
    $gene_level{$gene_level}="$gene_level{$gene_level}|$split1[0]";
    $st1=1;
    print out1 "$split1[0]\t$gene_level\t$arr1[0]";
    for (my $i=1; $i<=$#arr1; $i++){
     print out1 ":$arr1[$i]";
    }
    print out1 "\n";
    last LOOP1;
   }
  }

  if ($st1==0){
   my $temp=$#gene_level+1;
   $gene_level{$temp}=$split1[0];
   print out1 "$split1[0]\t$temp\t$arr1[0]";
   for (my $i=1; $i<=$#arr1; $i++){
    print out1 ":$arr1[$i]";
   }
   print out1 "\n";
  }
 }
}
close file1;
close out1;

open file1, "gunzip -c $sam |";
my (%pair_pos,%pair_flag,%pair_insert,%dedup_sr,%sr_pos,%samflag,%mq2,%cigar,%dedup_pe);
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 
 # paired-end analysis
 # MQ>=40, all insert sizes, mate on the same chromosome, both reads mapped, and primary alignment
 if ($split1[4]>=40 && $split1[6] eq "=" && !($split1[1]&4) && !($split1[1]&8) && !($split1[1]&256)){
  my $temp4="$split1[0].$split1[1].$split1[3]";

  my @split11=split /M/,$split1[5];
  my $ext=-1;
  foreach my $split11 (@split11){
   unless ($split11=~/[A-Za-z]$/){
    my @split12=split /[A-Za-z]/,$split11;
    $ext+=$split12[$#split12];
   }
  }
  my $read_end=$split1[3]+$ext;

  # paired-end reads opposite direction, forward read (read NOT reverse strand, mate reverse strand)
  if (!($split1[1]&16) && ($split1[1]&32)){
   $dedup_pe{$temp4}=1;
   $pair_pos{$split1[0]}="$split1[3]:$read_end|$split1[7]:$split1[7]";
   $pair_insert{$split1[0]}="$split1[8]";
   $samflag{$split1[0]}=$split1[1];
   $mq2{$split1[0]}=$split1[4];
   $cigar{$split1[0]}=$split1[5];

   my $flag1="FR_inward";
   if ($split1[8]<0){ #negative insert size
    $flag1="FR_outward";
   }
   $pair_flag{$split1[0]}=$flag1;
  }

  # paired-end reads both forward direction
  elsif (!($split1[1]&16) && !($split1[1]&32)){
   if ($split1[1]&64){ #only keep the first in pair
    $dedup_pe{$temp4}=1;
    $pair_pos{$split1[0]}="$split1[3]:$read_end|$split1[7]:$split1[7]";
    $pair_insert{$split1[0]}="$split1[8]";
    $pair_flag{$split1[0]}="FF";
    $samflag{$split1[0]}=$split1[1];
    $mq2{$split1[0]}=$split1[4];
    $cigar{$split1[0]}=$split1[5];
   }
  }

  # paired-end reads both reverse direction
  elsif (($split1[1]&16) && ($split1[1]&32)){
   if ($split1[1]&64){ #only keep the first in pair
    $dedup_pe{$temp4}=1;
    $pair_pos{$split1[0]}="$split1[3]:$read_end|$split1[7]:$split1[7]";
    $pair_insert{$split1[0]}="$split1[8]";
    $pair_flag{$split1[0]}="RR";
    $samflag{$split1[0]}=$split1[1];
    $mq2{$split1[0]}=$split1[4];
    $cigar{$split1[0]}=$split1[5];
   }
  }
 }

 # split-read analysis - supplementary alignment
 if ($split1[9]=~/^SA/ && $split1[4]>=40 && $split1[6] eq "="){
  my @split11=split /M/,$split1[5];
  my $ext=-1;
  foreach my $split11 (@split11){
   unless ($split11=~/[A-Za-z]$/){
    my @split12=split /[A-Za-z]/,$split11;
    $ext+=$split12[$#split12];
   }
  }
  my $read_end=$split1[3]+$ext;

  my @split2=split /\,/,$split1[9];
  my @split3=split /\:/,$split2[0];
  my $temp1="$split1[0].$split1[3]";
  unless (exists $dedup_sr{$temp1}){
   if ($chr eq $split3[$#split3]){
    my $temp2="$split1[0].$split2[1]";
    $dedup_sr{$temp2}=1;

    my @split13=split /M/,$split2[3];
    my $ext2=-1;
    foreach my $split13 (@split13){
     unless ($split13=~/[A-Za-z]$/){
      my @split14=split /[A-Za-z]/,$split13;
      $ext2+=$split14[$#split14];
     }
    }
    my $read_end2=$split2[1]+$ext2;
    
    my $temp3="$split1[0]_$split1[3]";
    $sr_pos{$temp3}="$split1[3]:$read_end|$split2[1]:$read_end2";
    $samflag{$temp3}=$split1[1];
    $mq2{$temp3}="$split1[4],$split2[4]";
    $cigar{$temp3}="$split1[5],$split2[3]";
   }
  }
 }
}

# update end position for paired-end reads (reverse strand or 2nd in pair)
open file1, "gunzip -c $sam |";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;

 # paired-end analysis
 # MQ>=40, all insert sizes, mate on the same chromosome, both reads mapped, and primary alignment
 if ($split1[4]>=40 && $split1[6] eq "=" && !($split1[1]&4) && !($split1[1]&8) && !($split1[1]&256)){
  my @split11=split /M/,$split1[5];
  my $ext=-1;
  foreach my $split11 (@split11){
   unless ($split11=~/[A-Za-z]$/){
    my @split12=split /[A-Za-z]/,$split11;
    $ext+=$split12[$#split12];
   }
  }
  my $read_end=$split1[3]+$ext;

  if (exists $pair_pos{$split1[0]}){
   my @split2=split /\|/,$pair_pos{$split1[0]};
   my @split3=split /\:/,$split2[1];
   my $temp5="$split1[0].$split1[1].$split1[3]";

   if ($split3[0] eq $split1[3]){
    unless (exists $dedup_pe{$temp5}){
     $dedup_pe{$temp5}=1;
     $pair_pos{$split1[0]}="$split2[0]|$split3[0]:$read_end";
     $samflag{$split1[0]}="$samflag{$split1[0]},$split1[1]";
     $mq2{$split1[0]}="$mq2{$split1[0]},$split1[4]";
     $cigar{$split1[0]}="$cigar{$split1[0]},$split1[5]";
    }
   }
  }
 }
}
close file1;

open out1, ">$dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script05_file2.txt";
my @pair_pos=keys %pair_pos;
@pair_pos=sort(@pair_pos);
foreach my $pair_pos (@pair_pos){
 my @split1=split /\|/,$pair_pos{$pair_pos};
 my @split2=split /\:/,$split1[0];
 my @split3=split /\:/,$split1[1];
 my @split4;
 if ($split2[0]<=$split3[0]){
  @split4=($split2[0],$split2[1],$split3[0],$split3[1]);
 }
 elsif ($split2[0]>$split3[0]){
  @split4=($split3[0],$split3[1],$split2[0],$split2[1]);
 }

 print out1 "$pair_pos\t$samflag{$pair_pos}\t$chr\t$split4[0]\t$split4[1]\t$split4[2]\t$split4[3]\t$pair_insert{$pair_pos}\t$mq2{$pair_pos}\t$cigar{$pair_pos}\t$pair_flag{$pair_pos}\n";
}

my @sr_pos=keys %sr_pos;
@sr_pos=sort(@sr_pos);
foreach my $sr_pos (@sr_pos){
 my @split1=split /\|/,$sr_pos{$sr_pos};
 my @split2=split /\:/,$split1[0];
 my @split3=split /\:/,$split1[1];
 my @split4=($split2[0],$split2[1],$split3[0],$split3[1]);
 @split4=sort{$a<=>$b}@split4;
 my $size=$split4[2]-$split4[1]+1;

 my @split5=split /\_/,$sr_pos;
 print out1 "$split5[0]\t$samflag{$sr_pos}\t$chr\t$split4[0]\t$split4[1]\t$split4[2]\t$split4[3]\t$size\t$mq2{$sr_pos}\t$cigar{$sr_pos}\tSR\n";
}
close out1;

system ("gzip -f $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script05_file1.txt $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script05_file2.txt $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script05_file3.txt");
exit 2;
