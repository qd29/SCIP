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

my $cwd=abs_path();
my ($dir,$cangene,$gnomad_common,$all_exons,$lowqual_reg,$lowqual_breakpoint);
open file1, "<$cwd/pipeline_config.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq "CURRENT_PATH"){
  $dir=$split1[1];
 }
 if ($split1[0] eq "gene_interest"){
  $cangene=$split1[1];
 }
 if ($split1[0] eq "gnomAD_common"){
  $gnomad_common=$split1[1];
 }
 if ($split1[0] eq "all_exons"){
  $all_exons=$split1[1];
 }
 if ($split1[0] eq "lowqual_reg"){
  $lowqual_reg=$split1[1];
 }
 if ($split1[0] eq "lowqual_bkp"){
  $lowqual_breakpoint=$split1[1];
 }
}
close file1;

my %can;
open file1, "<$cangene";
while (<file1>){
 chomp;
 $can{$_}=1;
}
close file1;

my @jct;
open file1, "<$all_exons";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq $chr){
  push @jct, $split1[1];
 }
}
close file1;

my $leftjct=0; my $rightjct=0;
foreach my $jct (@jct){
 if (abs($start-$jct)<=10){
  $leftjct=1;
 }
 if (abs($end-$jct)<=10){
  $rightjct=1;
 }
}

open out1, ">$dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script08_file1.txt";
open out2, ">$dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script08_file2.txt";
print out2 "HEADER\t0\n";
my @flag; my $clingen_reg_full=0;
if ($type eq "DEL" || $type eq "INS" || $type eq "BND"){
 my $clingen_hi_region=0; my $clingen_hi_gene=0; my $can_gene=0; my $pli9_loeuf35=0;
 my $pli5=0; my $dis_dom=0; my $dis_any=0;
 open file1, "gunzip -c $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script07_file1.txt.gz |";
 while (<file1>){
  chomp;
  my @split2=split /\t/,$_;
  my @split3=split /\>/,$split2[0];
  my @split4=split /\</,$split3[1];
  if (exists $can{$split4[0]}){
   if ($can_gene eq 0){
    $can_gene=$split4[0];
   }
   else{
    $can_gene="$can_gene, $split4[0]";
   }
  }
  if ($split2[7] ne "NA" && ($split2[7]>=0.9 || $split2[8]<=0.35)){
   if ($pli9_loeuf35 eq 0){
    $pli9_loeuf35=$split4[0];
   }
   else{
    $pli9_loeuf35="$pli9_loeuf35, $split4[0]";
   }
  }
  if ($split2[7] ne "NA" && $split2[7]>=0.5){
   if ($pli5 eq 0){
    $pli5=$split4[0];
   }
   else{
    $pli5="$pli5, $split4[0]";
   }
  }
  if ($split2[13] ne "NA" || $split2[14] ne "NA"){
   if ($split2[13] ne "NA"){
    my @split5=split /\<br\/\>/,$split2[13];
    foreach my $split5 (@split5){
     if ($dis_any eq 0){
      $dis_any=$split5;
     }
     else{
      $dis_any="$dis_any, $split5";
     }
    }
   }
   if ($split2[14] ne "NA"){
    my @split5=split /\<br\/\>/,$split2[14];
    foreach my $split5 (@split5){
     if ($dis_any eq 0){
      $dis_any=$split5;
     }
     else{
      $dis_any="$dis_any, $split5";
     }
    }
   }
  }
  if ($split2[17]==1){
   $dis_dom=1;
  }
 }
 close file1;

 open file1, "gunzip -c $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script06_file1.txt.gz |";
 while (<file1>){
  chomp;
  my @split2=split /\t/,$_;
  if ($split2[0] eq "ClinGen_Region" && ($split2[5]=~/Sufficient/ || $split2[5]=~/Emerging/ || $split2[5]=~/Little/)){
   # if the ClinGen item is a region, require it to be fully contained within the DEL
   if ($split2[1]=~/ISCA/ && $start<=$split2[3] && $end>=$split2[4]){
    if ($clingen_hi_region eq 0){
     $clingen_hi_region=$split2[1];
    }
    else{
     $clingen_hi_region="$clingen_hi_region, $split2[1]";
    }
   }
   unless ($split2[1]=~/ISCA/){
    if ($clingen_hi_gene eq 0){
     $clingen_hi_gene=$split2[1];
    }
    else{
     $clingen_hi_gene="$clingen_hi_gene, $split2[1]";
    }
   }
  }
 }
 close file1;

 my $order=99;
 if ($clingen_hi_region ne 0 || $clingen_hi_gene ne 0){
  $order=1;
 }
 elsif ($pli9_loeuf35 ne 0 || $can_gene ne 0){
  $order=2;
 }
 elsif ($pli5 ne 0){
  $order=3;
 }
 elsif ($dis_dom==1){
  $order=4;
 }
 elsif ($dis_any ne 0){
  $order=5;
 }

 print out1 "$proband.$chr.$start.$end.$type\t$order\t";
 if ($clingen_hi_region ne 0){
  push @flag, "ClinGen_HI_Region";
  print out2 "This $type is fully contained in ClinGen HI region(s) ($clingen_hi_region)\t1\n";
  $clingen_reg_full=1;
 }
 if ($clingen_hi_gene ne 0){
  push @flag, "ClinGen_HI_Gene";
  print out2 "This $type overlaps ClinGen HI gene(s) ($clingen_hi_gene)\t1\n";
 }
 if ($can_gene ne 0){
  push @flag, "Candidate_Gene";
  print out2 "This $type overlaps candidate gene(s) ($can_gene)\t1\n";
 }
 if ($pli9_loeuf35 ne 0){
  push @flag, "pLI_0.9_LOEUF_0.35";
  print out2 "This $type overlaps gene(s) with pLI >= 0.9 and/or LOEUF <= 0.35 ($pli9_loeuf35)\t1\n";
 }
 elsif ($pli5 ne 0){
  push @flag, "pLI_0.5";
  print out2 "This $type overlaps gene(s) with pLI >= 0.5 but < 0.9 ($pli5)\t1\n";
 }
 if ($dis_dom==1){
  push @flag, "Dom_Disease";
  print out2 "This $type overlaps gene(s) associated with dominant disorder(s): $dis_any\t1\n";
 }
 elsif ($dis_any ne 0){
  push @flag, "NonDom_Disease";
  print out2 "This $type overlaps gene(s) associated with only non-dominant disorder(s): $dis_any\t1\n";
 }
}

if ($type eq "DUP"){
 my $clingen_ts_full=0; my $clingen_ts_overlap=0; my $clingen_intragenic_hi=0; my $can_gene=0;
 my $intragenic_pli9_loeuf35=0; my $intragenic_pli5=0; my $intragenic_dis_dom=0; my $intragenic_dis_any=0;
 my $nonintragenic_dis_dom=0; my $nonintragenic_dis_any=0;
 open file1, "gunzip -c $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script07_file1.txt.gz |";
 my %intragenic;
 while (<file1>){
  chomp;
  my @split2=split /\t/,$_;
  my @split3=split /\>/,$split2[0];
  my @split4=split /\</,$split3[1];
  if (exists $can{$split4[0]}){
   if ($can_gene eq 0){
    $can_gene=$split4[0];
   }
   else{
    $can_gene="$can_gene, $split4[0]";
   }
  }
  if ($split2[4] eq "Intragenic"){
   $intragenic{$split4[0]}=1;
   if ($split2[7] ne "NA" && ($split2[7]>=0.9 || $split2[8]<=0.35)){
    if ($intragenic_pli9_loeuf35 eq 0){
     $intragenic_pli9_loeuf35=$split4[0];
    }
    else{
     $intragenic_pli9_loeuf35="$intragenic_pli9_loeuf35, $split4[0]";
    }
   }
   if ($split2[7] ne "NA" && $split2[7]>=0.5){
    if ($intragenic_pli5 eq 0){
     $intragenic_pli5=$split4[0];
    }
    else{
     $intragenic_pli5="$intragenic_pli5, $split4[0]";
    }
   }
  }

  if ($split2[13] ne "NA" || $split2[14] ne "NA"){
   if ($split2[4] eq "Intragenic"){
    if ($split2[13] ne "NA"){
     my @split5=split /\<br\/\>/,$split2[13];
     foreach my $split5 (@split5){
      if ($intragenic_dis_any eq 0){
       $intragenic_dis_any=$split5;
      }
      else{
       $intragenic_dis_any="$intragenic_dis_any, $split5";
      }
     }
    }
    if ($split2[14] ne "NA"){
     my @split5=split /\<br\/\>/,$split2[14];
     foreach my $split5 (@split5){
      if ($intragenic_dis_any eq 0){
       $intragenic_dis_any=$split5;
      }
      else{
       $intragenic_dis_any="$intragenic_dis_any, $split5";
      }
     }
    }
   }
   else{
    if ($split2[13] ne "NA"){
     my @split5=split /\<br\/\>/,$split2[13];
     foreach my $split5 (@split5){
      if ($nonintragenic_dis_any eq 0){
       $nonintragenic_dis_any=$split5;
      }
      else{
       $nonintragenic_dis_any="$nonintragenic_dis_any, $split5";
      }
     }
    }
    if ($split2[14] ne "NA"){
     my @split5=split /\<br\/\>/,$split2[14];
     foreach my $split5 (@split5){
      if ($nonintragenic_dis_any eq 0){
       $nonintragenic_dis_any=$split5;
      }
      else{
       $nonintragenic_dis_any="$nonintragenic_dis_any, $split5";
      }
     }
    }
   }
  }
  if ($split2[17]==1){
   if ($split2[4] eq "Intragenic"){
    $intragenic_dis_dom=1;
   }
   else{
    $nonintragenic_dis_dom=1;
   }
  }
 } 
 close file1;

 open file1, "gunzip -c $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script06_file1.txt.gz |";
 while (<file1>){
  chomp;
  my @split2=split /\t/,$_;
  if ($split2[0] eq "ClinGen_Region"){
   my @split3=split /\>/,$split2[1];
   my @split4=split /\</,$split3[1];
   if ($split2[6]=~/Sufficient/ || $split2[6]=~/Emerging/ || $split2[6]=~/Little/){
    if ($start<=$split2[3] && $end>=$split2[4]){
     if ($clingen_ts_full eq 0){
      $clingen_ts_full=$split2[1];
     }
     else{
      $clingen_ts_full="$clingen_ts_full, $split2[1]";
     }
    }
    unless ($split2[1]=~/ISCA/){
     #genes only
     if ($clingen_ts_overlap eq 0){
      $clingen_ts_overlap=$split2[1];
     }
     else{
      $clingen_ts_overlap="$clingen_ts_overlap, $split2[1]";
     }
    }
   }

   if ($split2[5]=~/Sufficient/ || $split2[5]=~/Emerging/ || $split2[5]=~/Little/){
    if (exists $intragenic{$split4[0]}){
     if ($clingen_intragenic_hi eq 0){
      $clingen_intragenic_hi=$split2[1];
     }
     else{
      $clingen_intragenic_hi="$clingen_intragenic_hi, $split2[1]";
     }
    }
   }
  }
 }
 close file1;

 my $order=99;
 if ($clingen_ts_full ne 0 || $clingen_ts_overlap ne 0 || $clingen_intragenic_hi ne 0){
  $order=1;
 }
 elsif ($intragenic_pli9_loeuf35 ne 0 || $can_gene ne 0){
  $order=2;
 }
 elsif ($intragenic_pli5 ne 0){
  $order=3;
 }
 elsif ($intragenic_dis_dom==1){
  $order=4;
 }
 elsif ($nonintragenic_dis_dom==1){
  $order=5;
 }
 elsif ($intragenic_dis_any ne 0){
  $order=6;
 }
 elsif ($nonintragenic_dis_any ne 0){
  $order=7;
 }
 print out1 "$proband.$chr.$start.$end.$type\t$order\t";
 if ($clingen_ts_full ne 0){
  push @flag, "ClinGen_TS_Contained";
  print out2 "This DUP fully contains ClinGen TS region(s) or gene(s) ($clingen_ts_full)\t1\n";
  $clingen_reg_full=1;
 } 
 if ($clingen_ts_overlap ne 0){
  push @flag, "ClinGen_TS_Gene_Overlap";
  print out2 "This DUP overlaps ClinGen TS gene(s) ($clingen_ts_overlap)\t1\n";
 }
 if ($clingen_intragenic_hi ne 0){
  push @flag, "ClinGen_Intragenic_HI_Gene";
  print out2 "This intragenic DUP overlaps ClinGen HI gene(s) ($clingen_intragenic_hi)\t1\n";
 }
 if ($can_gene ne 0){
  push @flag, "Candidate_Gene";
  print out2 "This DUP overlaps candidate gene(s) ($can_gene)\t1\n";
 }
 if ($intragenic_pli9_loeuf35 ne 0){
  push @flag, "Intragenic_pLI_0.9_LOEUF_0.35";
  print out2 "This intragenic DUP overlaps gene(s) with pLI >= 0.9 and/or LOEUF <= 0.35 ($intragenic_pli9_loeuf35)\t1\n";
 }
 elsif ($intragenic_pli5 ne 0){
  push @flag, "Intragenic_pLI_0.5";
  print out2 "This intragenic DUP overlaps gene(s) with pLI >= 0.5 ($intragenic_pli5)\t1\n";
 }
 if ($intragenic_dis_dom==1){
  push @flag, "Intragenic_Dom_Disease";
  print out2 "This intragenic DUP overlaps gene(s) associated with dominant disorder(s): $intragenic_dis_any\t1\n";
 }
 elsif ($intragenic_dis_any ne 0){
  push @flag, "Intragenic_NonDom_Disease";
  print out2 "This intragenic DUP overlaps gene(s) associated with only non-dominant disorder(s): $intragenic_dis_any\t1\n";
 }
 if ($nonintragenic_dis_dom==1){
  push @flag, "NonIntragenic_Dom_Disease";
  print out2 "This non-intragenic DUP overlaps gene(s) associated with dominant disorder(s): $nonintragenic_dis_any\t1\n";
 }
 elsif ($nonintragenic_dis_any ne 0){
  push @flag, "NonIntragenic_NonDom_Disease";
  print out2 "This non-intragenic DUP overlaps gene(s) associated with only non-dominant disorder(s): $nonintragenic_dis_any\t1\n";
 }
}

my $neg_info=0;
open file1, "<$gnomad_common";
my %gnomad_dedup;
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($chr eq $split1[0] && $type eq $split1[3] && $split1[2]>=$start && $split1[1]<=$end){
  for (my $i=$split1[1]; $i<=$split1[2]; $i++){
   if ($i>=$start && $i<=$end){
    $gnomad_dedup{$i}=1;
   }
  }
 } 
}
close file1;
my @gnomad_dedup=keys %gnomad_dedup;
my $gnomad_overlap=($#gnomad_dedup+1)/($end-$start+1);
if ($gnomad_overlap>=0.5){
 $neg_info=1;
 my $gnomad_overlap_percent=int($gnomad_overlap*1e2*1e2)/1e2;
 print out2 "This $type has $gnomad_overlap_percent% overlap with gnomAD $type(s) with popmax >= 1%\t2\n";
 push @flag, "HIGH_gnomAD_SV_common_var_OVERLAP";
 if ($gnomad_overlap==1){
  push @flag, "FULL_gnomAD_SV_common_var_OVERLAP";
 }
}

if ($leftjct==1 && $rightjct==1 && $type eq "DEL"){
 $neg_info=1;
 push @flag, "BOTH_BND_NEAR_JCT_DEL";
 print out2 "This $type has both breakpoints within 10 bp of an intron-exon junction - rule out pseudogene insertion\t2\n";
}

my %erdsp_dedup; my $manta_break=0; my $manta_ext=5000;
open file1, "<$lowqual_reg";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[3] eq $type || $split1[3] eq "ALL" || $split1[3] eq "ALL_SR"){
  if ($chr eq $split1[0] && $split1[2]>=$start && $split1[1]<=$end){
   for (my $i=$split1[1]; $i<=$split1[2]; $i++){
    if ($i>=$start && $i<=$end){
     $erdsp_dedup{$i}=1;
    }
   }
  }
 }
}
close file1;

open file1, "<$lowqual_breakpoint";
LOOP1: while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($chr eq $split1[0] && $type eq $split1[5] && $start>=$split1[1]-$manta_ext && $start<=$split1[2]+$manta_ext && $end>=$split1[3]-$manta_ext && $end<=$split1[4]+$manta_ext){
  $manta_break=1; last LOOP1;
 }
}
close file1;

my @erdsp_dedup=keys %erdsp_dedup;
my $erdsp_overlap=($#erdsp_dedup+1)/($end-$start+1);
if ($erdsp_overlap>=0.5){
 $neg_info=1;
 my $erdsp_overlap_percent=int($erdsp_overlap*1e2*1e2)/1e2;
 print out2 "This $type has $erdsp_overlap_percent% overlap with ERDS+ $type recurrence region(s) or other low quality region(s)\t2\n";
 push @flag, "HIGH_ERDS+_OVERLAP";
}
if ($manta_break==1){
 $neg_info=1;
 print out2 "This $type has both breakpoints within $manta_ext bp of Manta $type recurrence region(s)\t2\n";
 push @flag, "HIGH_MANTA_OVERLAP";
}

my @pext;
open file1, "gunzip -c $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script05_file1.txt.gz |";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[3] ne "NA" && $split1[0]>=$start && $split1[0]<=$end){
  push @pext, $split1[3];
 }
}
close file1;
@pext=sort{$a<=>$b}@pext;
my $med_pext=0; my $max_pext=0;
if ($#pext>=0){
 $med_pext=$pext[int($#pext/2)];
 $max_pext=$pext[$#pext];
}
#if ($med_pext<=0.4){
if ($max_pext<=0.4){
 #$neg_info=1;
 #my $med_pext_round=int($med_pext*1e2)/1e2;
 my $max_pext_round=int($max_pext*1e2)/1e2;
 print out2 "The maximum pext within the $type is low ($max_pext_round) - double check: gene may have a highly expressed non-coding transcript or not expressed in GTEx tissues\t2\n";
 push @flag, "LOW_PEXT";
}
if ($max_pext<=0.2){
 push @flag, "MAX_PEXT_LESS_THAN_0.2";
}

# full in ClinGen regions overrides neg info
if ($clingen_reg_full==1){
 $neg_info=0;
}
print out1 "$neg_info\t";

my $anomalous=0;
open file1, "$dir/d1stat/$proband.$chr.$start.$end.$type.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[6]>0){
  $anomalous=1;
  push @flag, "PAIR_READ_SUPPORT";
  print out2 "This $type is supported by $split1[6] anomalous read pair(s)\t1\n";
 } 
 if ($split1[9]>0){
  $anomalous=1;
  push @flag, "SPLIT_READ_SUPPORT";
  print out2 "This $type is supported by $split1[9] anomalous split-read(s)\t1\n";
 }
 if ($split1[6]==0 && $split1[9]==0 && ($type eq "DEL" || $type eq "DUP")){
  print out2 "There are no anomalous read(s) supporting this $type\t2\n";
 }
 last;
}
close file1;
print out1 "$anomalous\t";

if ($#flag>=0){
 print out1 "$flag[0]";
 for (my $i=1; $i<=$#flag; $i++){
  print out1 " $flag[$i]";
 }
}
print out1 "\n";
close out1;
close out2;

system ("gzip -f $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script08_file1.txt $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script08_file2.txt");
exit 2;
