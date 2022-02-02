#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Cwd 'abs_path';
my %opts;
getopt ('c:s:e:p:t:',\%opts);
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

# default window size is 5 kb, however, if the CNV is < 5 kb in length, set window size to 100 bp
my $winsize=5000;
if ($exp_length<=5000){
 $winsize=100;
}

my $cwd=abs_path();
my ($sampleid_file,$dir,$prefix,$aln_dir);
open file1, "<$cwd/pipeline_config_hg38.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq "SAMPLE_ID"){
  $sampleid_file=$split1[1];
 }
 if ($split1[0] eq "CURRENT_PATH"){
  $dir=$split1[1];
 }
 if ($split1[0] eq "PREFIX"){
  $prefix=$split1[1];
 }
 if ($split1[0] eq "ALIGNMENT_PATH"){
  $aln_dir=$split1[1];
 }
}
close file1;
my @prefix=split /\;/,$prefix;

my @cyc=("d1temp_server");

open file1, "<$sampleid_file";
my %id;
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 my @split2=split /\-/,$split1[0];
 $id{$split1[1]}=$split1[0];
}
close file1;

my $preload=0;
my ($depth,$sam);
foreach my $cyc (@cyc){
 if (-e "$dir/$cyc/$proband.$chr.$start_original.$end_original.$type.SAM.txt.gz"){
  $preload=1;
  $depth="$dir/$cyc/$proband.$chr.$start_original.$end_original.$type.DEPTH.txt.gz";
  $sam="$dir/$cyc/$proband.$chr.$start_original.$end_original.$type.SAM.txt.gz";
 }
}

my $file;
if ($preload==0){
 my $st=0;
 if (-e "$aln_dir/$id{$proband}/$id{$proband}.cram"){
  $file="$aln_dir/$id{$proband}/$id{$proband}.cram"; $st=1;
 }
 elsif (-e "$aln_dir/$id{$proband}/$id{$proband}.bam"){
  $file="$aln_dir/$id{$proband}/$id{$proband}.bam"; $st=1;
 }
 foreach my $prefix (@prefix){
  if (-e "$aln_dir/$id{$proband}/$id{$proband}.$prefix.cram"){
   $file="$aln_dir/$id{$proband}/$id{$proband}.$prefix.cram"; $st=1;
  }
  elsif (-e "$aln_dir/$id{$proband}/$id{$proband}.$prefix.bam"){
   $file="$aln_dir/$id{$proband}/$id{$proband}.$prefix.bam"; $st=1;
  }
 }
 if ($st==0){
  die "No BAM/CRAM file found for $id{$proband}!";
 }

 print "Generating new temporary SAM/DEPTH files\n";
 system ("samtools view -F 0x400 \"$file\" chr$chr:$start-$end |cut -f1,2,3,4,5,6,7,8,9,12 -d\$'\t' > $dir/d1temp_server/$proband.$chr.$start_original.$end_original.$type.SAM.txt");
 system ("samtools depth \"$file\" -r chr$chr:$start-$end > $dir/d1temp_server/$proband.$chr.$start_original.$end_original.$type.DEPTH.txt");
 system ("gzip -f $dir/d1temp_server/$proband.$chr.$start_original.$end_original.$type.SAM.txt $dir/d1temp_server/$proband.$chr.$start_original.$end_original.$type.DEPTH.txt");
 $depth="$dir/d1temp_server/$proband.$chr.$start_original.$end_original.$type.DEPTH.txt.gz";
 $sam="$dir/d1temp_server/$proband.$chr.$start_original.$end_original.$type.SAM.txt.gz";
}
else{
 print "Reusing existing temporary SAM/DEPTH files\n";
}

my %depth;
open file1, "gunzip -c $depth |";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 $depth{$split1[1]}=$split1[2];
}
close file1;

my (%gq,%pair,%split,%largeinsert,%splitread,%dedup_sr);
open file1, "gunzip -c $sam |";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if (exists $gq{$split1[3]}){
  $gq{$split1[3]}="$gq{$split1[3]}|$split1[4]";
 }
 else{
  $gq{$split1[3]}=$split1[4];
 }

 # paired-end analysis
 # MQ>=40, insert size is 80% to 120% of expected, mate on the same chromosome, and both reads mapped
 if ($split1[4]>=40 && abs($split1[8])>=0.8*$exp_length && abs($split1[8])<=1.2*$exp_length && $split1[6] eq "=" && !($split1[1]&4) && !($split1[1]&8)){
  # paired-end reads opposite direction, forward read (read NOT reverse strand, mate reverse strand)
  my $st=0;
  if (!($split1[1]&16) && ($split1[1]&32)){
   $st=1;
   my $flag1="FR_inward";
   if ($split1[8]<0){ #negative insert size
    $flag1="FR_outward";
   }

   if (exists $pair{$split1[3]}){
    $pair{$split1[3]}="$pair{$split1[3]}|$split1[7]:$flag1";
   }
   else{
    $pair{$split1[3]}="$split1[7]:$flag1";
   }
  }

  # paired-end reads both forward direction
  elsif (!($split1[1]&16) && !($split1[1]&32)){
   if ($split1[1]&64){ #only keep the first in pair
    $st=1;
    if (exists $pair{$split1[3]}){
     $pair{$split1[3]}="$pair{$split1[3]}|$split1[7]:FF";
    }
    else{
     $pair{$split1[3]}="$split1[7]:FF";
    }
   }
  }

  # paired-end reads both reverse direction
  elsif (($split1[1]&16) && ($split1[1]&32)){
   if ($split1[1]&64){ #only keep the first in pair
    $st=1;
    if (exists $pair{$split1[3]}){
     $pair{$split1[3]}="$pair{$split1[3]}|$split1[7]:RR";
    }
    else{
     $pair{$split1[3]}="$split1[7]:RR";
    }
   }
  }

  # paired-end reads opposite direction, reverse read, no need to output
  elsif (($split1[1]&16) && !($split1[1]&32)){}
  else{
   print "WARNING $proband $start_original $end_original $type unrecognized flag $split1[1]\n";
  }

  if ($st==1){
   $largeinsert{$split1[3]}++;
   $largeinsert{$split1[7]}++;
  }
 }

 # split-read analysis - supplementary alignment
 if ($split1[9]=~/^SA/ && $split1[4]>=40 && $split1[6] eq "="){
  my @split2=split /\,/,$split1[9];
  my @split3=split /\:/,$split2[0];
  my $temp1="$split1[0].$split1[3]";
  unless (exists $dedup_sr{$temp1}){
   my $st2=0;
   if (abs($split2[1]-$split1[3])>=0.8*$exp_length && abs($split2[1]-$split1[3])<=1.2*$exp_length){
    $st2=1;
   }
   if (abs(abs($split2[1]-$split1[3])-$exp_length)<=200){
    $st2=1;
   }

   if ($chr eq $split3[$#split3] && $st2==1){
    my $temp2="$split1[0].$split2[1]";
    $dedup_sr{$temp2}=1;
    if (exists $split{$split1[3]}){
     $split{$split1[3]}="$split{$split1[3]}|$split2[1]";
    }
    else{
     $split{$split1[3]}=$split2[1];
    }
    $splitread{$split1[3]}++;
    $splitread{$split2[1]}++;
   }
  }
 } 
}
close file1;

open out1, ">$dir/$proband.$chr.$start.$end.temp1.txt";
for (my $i=$start; $i<=$end; $i+=$winsize){
 my $ct_depth=0; my $tot_depth=0;
 my $ct_gq=0; my $tot_gq=0;
 my $tot_largeinsert=0; my $tot_splitread=0;
 my $right=$i+$winsize;

 for (my $j=$i; $j<$i+$winsize; $j++){
  if (exists $depth{$j}){
   $ct_depth++; $tot_depth+=$depth{$j};
  }
  if (exists $gq{$j}){
   my @split1=split /\|/,$gq{$j};
   foreach my $split1 (@split1){
    $ct_gq++; $tot_gq+=$split1;
   }
  }
  if (exists $largeinsert{$j}){
   $tot_largeinsert+=$largeinsert{$j};
  }
  if (exists $splitread{$j}){
   $tot_splitread+=$splitread{$j};
  }
 }

 my $rt_depth="NA"; my $rt_gq="NA";
 if ($ct_depth>0){ 
  $rt_depth=$tot_depth/$ct_depth;
 }
 if ($ct_gq>0){
  $rt_gq=$tot_gq/$ct_gq;
 }
 my $flag=0;
 if ($i>=$start_original && $i<=$end_original){
  $flag=1;
 }
 print out1 "$chr\t$i\t$right\t$rt_depth\t$rt_gq\t$tot_largeinsert\t$tot_splitread\t$flag\n";
}
close out1;

open out1, ">$dir/$proband.$chr.$start.$end.temp2.txt";
my $inpair=0;
my @pair=keys %pair;
@pair=sort{$a<=>$b}@pair;
foreach my $pair (@pair){
 $inpair=1;
 my @split1=split /\|/,$pair{$pair};
 foreach my $split1 (@split1){
  my @split2=split /\:/,$split1;
  if ($split2[0]>$pair){
   print out1 "$pair\t$split2[0]\t$split2[1]\n";
  }
  else{
   print out1 "$split2[0]\t$pair\t$split2[1]\n";
  }
 }
}
close out1;

open out1, ">$dir/$proband.$chr.$start.$end.temp3.txt";
my $insplit=0;
my @split=keys %split;
@split=sort{$a<=>$b}@split;
foreach my $split (@split){
 $insplit=1;
 my @split1=split /\|/,$split{$split};
 foreach my $split1 (@split1){
  if ($split>$split1){
   print out1 "$split1\t$split\n";
  }
  else{
   print out1 "$split\t$split1\n";
  }
 }
}
close out1;

open out1, ">$dir/$proband.$chr.$start.$end.temp4.txt";
print out1 "x=read.table (\"$dir/$proband.$chr.$start.$end.temp1.txt\")\n";
if ($inpair==1){
 print out1 "y=read.table (\"$dir/$proband.$chr.$start.$end.temp2.txt\")\n";
}
else{
 print out1 "y=c()\n";
}
if ($insplit==1){
 print out1 "z=read.table (\"$dir/$proband.$chr.$start.$end.temp3.txt\")\n";
}
else{
 print out1 "z=c()\n";
}
print out1 "proband=\"$proband\"\nchr=\"$chr\"\nlb=$start_original\nrb=$end_original\ntype=\"$type\"\n";
open file1, "<$dir/SCIP_pri_03.R";
while (<file1>){
 print out1 $_;
}
close file1;
print out1 "write.table (cbind(median_depth_in,median_depth_out,median_depth_in/median_depth_out,median_mq_in,median_mq_out,median_mq_in/median_mq_out,ct_paired_end_support,ct_paired_end_opposite,ct_paired_end_non_supporting,ct_split_read),file=\"$dir/d1stat/$proband.$chr.$start_original.$end_original.$type.txt\",sep=\"\t\",col.names=F,row.names=F)\n";

close out1;
system ("R CMD BATCH --vanilla $dir/$proband.$chr.$start.$end.temp4.txt");
#system ("rm $dir/$proband.$chr.$start.$end.temp1.txt $dir/$proband.$chr.$start.$end.temp2.txt $dir/$proband.$chr.$start.$end.temp3.txt $dir/$proband.$chr.$start.$end.temp4.txt $dir/$proband.$chr.$start.$end.temp4.txt.Rout");
system ("rm $dir/$proband.$chr.$start.$end.temp1.txt $dir/$proband.$chr.$start.$end.temp2.txt $dir/$proband.$chr.$start.$end.temp3.txt $dir/$proband.$chr.$start.$end.temp4.txt");
system ("rm -f $dir/$proband.$chr.$start.$end.temp4.txt.Rout ~/$proband.$chr.$start.$end.temp4.txt.Rout");

exit 2;
