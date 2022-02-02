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
my ($clinvar_cnv,$gnomad_sv,$clingen_region,$cgc_cnv,$clingen_gene,$dir);
open file1, "<$cwd/pipeline_config.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq "ClinVar_CNV"){
  $clinvar_cnv=$split1[1];
 }
 if ($split1[0] eq "gnomAD_SV"){
  $gnomad_sv=$split1[1];
 }
 if ($split1[0] eq "ClinGen_dosage_region"){
  $clingen_region=$split1[1];
 }
 if ($split1[0] eq "ClinGen_dosage_gene"){
  $clingen_gene=$split1[1];
 }
 if ($split1[0] eq "cohort_CNV"){
  $cgc_cnv=$split1[1];
 }
 if ($split1[0] eq "CURRENT_PATH"){
  $dir=$split1[1];
 }
}
close file1;

open out1, ">$dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script06_file1.txt";
print out1 "HEADER\n";
open file1, "<$clinvar_cnv";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq $chr && $split1[2]>=$start && $split1[1]<=$end){
  my $in=1;
  if (($type eq "DEL" && $split1[4] eq "DUP") || ($type eq "DUP" && $split1[4] eq "DEL")){
   $in=0;
  }
  if ($in==1){
   print out1 "ClinVar\t$split1[0]\t$split1[1]\t$split1[2]\t$split1[3]\t$split1[4]\t$split1[5]\t$split1[6]\t\"$split1[8]\"\t";
   my @split2=split /\,/,$split1[7];
   foreach my $split2 (@split2){
    print out1 "<a href='https://www.ncbi.nlm.nih.gov/clinvar/$split2' target='_blank'>$split2</a> ";
   }
   my $path="OTHER";
   if ($split1[3] eq "Pathogenic" || $split1[3] eq "Likely pathogenic" || $split1[3] eq "Pathogenic/Likely pathogenic"){
    $path="P";
   }
   if ($split1[3] eq "Benign" || $split1[3] eq "Likely benign" || $split1[3] eq "Benign/Likely benign"){
    $path="B";
   }
   print out1 "\t$path\n";
  }
 }
}
close file1;

system ("tabix $gnomad_sv $chr:$start-$end > $dir/$proband.$chr.$start.$end.$type.d102temp1");
open file1, "<$dir/$proband.$chr.$start.$end.$type.d102temp1";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 my @split2=split /gnomAD\-SV\_v2\.1\_/,$split1[2];
 my @split3=split /END\=/,$split1[7]; my @split31=split /\;/,$split3[1];
 my @split4=split /SVTYPE\=/,$split1[7]; my @split41=split /\;/,$split4[1];
 my $popmax="NA";
 my @split5=split /POPMAX_AF\=/,$split1[7];
 if ($split5[1]){
  my @split51=split /\;/,$split5[1];
  $popmax=$split51[0];
 }
 my $in=1;
 if (($type eq "DEL" && $split41[0] eq "DUP") || ($type eq "DUP" && $split41[0] eq "DEL")){
  $in=0;
 }
 if ($type ne "INS" && $split41[0] eq "INS"){
  $in=0;
 }

 unless ($split1[6]=~"UNRESOLVED"){
  if ($in==1){
   print out1 "gnomAD_SV\t<a href='https://gnomad.broadinstitute.org/variant/$split2[1]?dataset=gnomad_sv_r2_1' target='_blank'>$split2[1]</a>\t$split1[0]\t$split1[1]\t$split31[0]\t$split41[0]\t$popmax\t$split1[6]\t$split2[1]\n";
  }
 }
}
close file1;
system ("rm $dir/$proband.$chr.$start.$end.$type.d102temp1");

open file1, "<$clingen_region";
my %order_clingen;
while (<file1>){
 chomp;
 if ($_=~/^#ISCA/){
  my @split1=split /\t/,$_;
  for (my $i=0; $i<=$#split1; $i++){
   $order_clingen{$split1[$i]}=$i;
  }
  last;
 }
}

my %clingen_score;
$clingen_score{"3"}="Sufficient (3)"; $clingen_score{"2"}="Emerging (2)";
$clingen_score{"1"}="Little (1)"; $clingen_score{"0"}="No evidence (0)";
$clingen_score{"30"}="Recessive (30)"; $clingen_score{"40"}="Unlikely (40)"; 

while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[$order_clingen{"Genomic Location"}] ne "tbd"){
  my @split2=split /\:/,$split1[$order_clingen{"Genomic Location"}];
  my @split3=split /\-/,$split2[1];
  my @split4=split /chr/,$split2[0];
  if ($split4[1] eq $chr && $split3[1]>=$start && $split3[0]<=$end){
   print out1 "ClinGen_Region\t<a href='https://dosage.clinicalgenome.org/clingen_region.cgi?id=$split1[0]' target='_blank'>$split1[0]</a>";
   my $hi_score="NA"; my $ts_score="NA";
   if (exists $clingen_score{$split1[$order_clingen{"Haploinsufficiency Score"}]}){$hi_score=$clingen_score{$split1[$order_clingen{"Haploinsufficiency Score"}]}}
   if (exists $clingen_score{$split1[$order_clingen{"Triplosensitivity Score"}]}){$ts_score=$clingen_score{$split1[$order_clingen{"Triplosensitivity Score"}]}}
   print out1 "\t$split4[1]\t$split3[0]\t$split3[1]\t$hi_score\t$ts_score\n";
  }
 }
}
close file1;

open file1, "<$clingen_gene";
my %order_clingen_gene;
while (<file1>){
 chomp;
 if ($_=~/^#Gene/){
  my @split1=split /\t/,$_;
  for (my $i=0; $i<=$#split1; $i++){
   $order_clingen_gene{$split1[$i]}=$i;
  }
  last;
 }
}

while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 my @split2=split /\:/,$split1[$order_clingen_gene{"Genomic Location"}];
 my @split3=split /\-/,$split2[1];
 my @split4=split /chr/,$split2[0];
 if ($split4[1] eq $chr && $split3[1]>=$start && $split3[0]<=$end){
  print out1 "ClinGen_Region\t<a href='https://dosage.clinicalgenome.org/clingen_gene.cgi?sym=$split1[0]' target='_blank'>$split1[0]</a>";
  my $hi_score="NA"; my $ts_score="NA";
  if (exists $clingen_score{$split1[$order_clingen_gene{"Haploinsufficiency Score"}]}){$hi_score=$clingen_score{$split1[$order_clingen_gene{"Haploinsufficiency Score"}]}}
  if (exists $clingen_score{$split1[$order_clingen_gene{"Triplosensitivity Score"}]}){$ts_score=$clingen_score{$split1[$order_clingen_gene{"Triplosensitivity Score"}]}}
  print out1 "\t$split4[1]\t$split3[0]\t$split3[1]\t$hi_score\t$ts_score\n";
 }
}
close file1;

system ("grep ^$chr $cgc_cnv |grep $type > $dir/$proband.$chr.$start.$end.$type.d102temp2");
open file1, "<$dir/$proband.$chr.$start.$end.$type.d102temp2";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[3] eq $type && $split1[4] ne $proband && $split1[0] eq $chr && $split1[2]>=$start && $split1[1]<=$end){
  my @split2=split /\-/,$split1[4];
  my @split3=split /\-/,$proband;
  my $temp1="$split2[0]-$split2[1]"; my $temp2="$split3[0]-$split3[1]";
  my $same_family_cgc=0;
  if ($temp1 eq $temp2){
   $same_family_cgc=1;
  }

  print out1 "cohort_CNV\t$_\t$same_family_cgc\n";
 }
}
close file1;
close out1;
system ("rm $dir/$proband.$chr.$start.$end.$type.d102temp2");
system ("gzip -f $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script06_file1.txt");

exit 2;
