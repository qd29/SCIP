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
my ($gnomad_constraints,$expression,$goterms,$gencc,$omim,$clingen_dosage_gene,$coding_exons,$dir,$strand,$search,$transcript_info,$noncoding);
open file1, "<$cwd/pipeline_config.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[0] eq "ClinGen_dosage_gene"){
  $clingen_dosage_gene=$split1[1];
 }
 if ($split1[0] eq "gnomAD_constraints"){
  $gnomad_constraints=$split1[1];
 }
 if ($split1[0] eq "expression_file"){
  $expression=$split1[1];
 }
 if ($split1[0] eq "GO_terms"){
  $goterms=$split1[1];
 }
 if ($split1[0] eq "OMIM"){
  $omim=$split1[1];
 }
 if ($split1[0] eq "GenCC"){
  $gencc=$split1[1];
 }
 if ($split1[0] eq "coding_exons"){
  $coding_exons=$split1[1];
 }
 if ($split1[0] eq "CURRENT_PATH"){
  $dir=$split1[1];
 }
 if ($split1[0] eq "gene_strand"){
  $strand=$split1[1];
 }
 if ($split1[0] eq "search_terms"){
  $search=$split1[1];
 }
 if ($split1[0] eq "transcript_info"){
  $transcript_info=$split1[1];
 }
 if ($split1[0] eq "noncoding"){
  $noncoding=$split1[1];
 }
} 
close file1;

my %noncoding_exc;
open file1, "<$noncoding";
while (<file1>){
 chomp;
 $noncoding_exc{$_}=1;
}
close file1;

my (%enst,%exon_ct,%coding_exon,%all_exon,%mane);
open file1, "<$transcript_info";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 $enst{$split1[0]}=$split1[1];
 $exon_ct{$split1[0]}=$split1[2];
 $coding_exon{$split1[0]}=$split1[3];
 $all_exon{$split1[0]}=$split1[4];
 $mane{$split1[0]}=$split1[5]; 
}
close file1;

my %strand;
open file1, "<$strand";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[1] eq 1){
  $strand{$split1[0]}="Forward";
 }
 elsif ($split1[1] eq -1){
  $strand{$split1[0]}="Reverse";
 }
}
close file1;

my (%hi,%ts,%clingen_score);
$clingen_score{'40'}="Unlikely"; $clingen_score{'30'}="Recessive";
$clingen_score{'3'}="Sufficient"; $clingen_score{'2'}="Emerging";
$clingen_score{'1'}="Little";
open file1, "<$clingen_dosage_gene";
my %order1;
while (<file1>){
 chomp;
 if ($_=~/^#Gene/){
  my @split1=split /\t/,$_;
  for (my $i=0; $i<=$#split1; $i++){
   $order1{$split1[$i]}=$i;
  }
  last;
 }
}

while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[$order1{'Haploinsufficiency Score'}]>0){
  $hi{$split1[0]}=$clingen_score{$split1[$order1{'Haploinsufficiency Score'}]};
 }
 if ($split1[$order1{'Triplosensitivity Score'}] ne "Not yet evaluated" && $split1[$order1{'Triplosensitivity Score'}]>0){
  $ts{$split1[0]}=$clingen_score{$split1[$order1{'Triplosensitivity Score'}]};
 }
}
close file1;

open file1, "gunzip -c $gnomad_constraints |";
my $dump2=<file1>; chomp $dump2; my %order_gnomad;
my @split4=split /\t/,$dump2;
for (my $i=0; $i<=$#split4; $i++){
 $order_gnomad{$split4[$i]}=$i;
}

my (%pli,%prec,%pnull,%loeuf);
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[$order_gnomad{'canonical'}] eq "true"){
  $pli{$split1[$order_gnomad{'gene_id'}]}=$split1[$order_gnomad{'pLI'}];
  $prec{$split1[$order_gnomad{'gene_id'}]}=$split1[$order_gnomad{'pRec'}];
  $pnull{$split1[$order_gnomad{'gene_id'}]}=$split1[$order_gnomad{'pNull'}];
  $loeuf{$split1[$order_gnomad{'gene_id'}]}=$split1[$order_gnomad{"oe_lof_upper"}];
 }
}
close file1;

open file1, "<$omim";
my %moi;
$moi{"Autosomal dominant"}="AD"; $moi{"Autosomal recessive"}="AR";
$moi{"X-linked dominant"}="XLD"; $moi{"X-linked recessive"}="XLR";
$moi{"X-linked"}="XL";

my (%omim,%domgene,%omimgene);
my $dump4=<file1>; $dump4=<file1>; $dump4=<file1>; $dump4=<file1>;
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[10] && $split1[12] && $split1[10] ne "" && $split1[12] ne ""){
  if ($split1[5] ne ""){
   $omimgene{$split1[10]}=$split1[5];
  }

  my @split2=split /\;\s/,$split1[12];
  foreach my $split2 (@split2){
   my @split3=split /\,\s/,$split2;
   my $mimid="NA";
   for (my $i=0; $i<=$#split3; $i++){
    if ($split3[$i]=~/^[0-9][0-9][0-9][0-9][0-9][0-9]/){
     $mimid=$i;
    }
   }
   if ($mimid ne "NA"){
    my $name=$split3[0];
    for (my $i=1; $i<$mimid; $i++){
     $name="$name, $split3[$i]";
    }
    my @split4=split /\s/,$split3[$mimid];
    my $moi="unknown";
    if ($#split3>$mimid){
     if (exists $moi{$split3[$mimid+1]}){
      $moi=$moi{$split3[$mimid+1]};
     }
     else{
      $moi=$split3[$mimid+1];
     }

     for (my $i=$mimid+2; $i<=$#split3; $i++){
      if (exists $moi{$split3[$i]}){
       $moi="$moi,$moi{$split3[$i]}";
      }
      else{
       $moi="$moi,$split3[$i]";
      }
     }
    }

    if ($moi=~/AD/ || $moi=~/XLD/ || $moi=~/dominant/){
     $domgene{$split1[10]}=1;
    }

    my $temp="$name ($moi)|$split4[0]";
    if (exists $omim{$split1[10]}){
     $omim{$split1[10]}="$omim{$split1[10]}::$temp"; 
    } 
    else{
     $omim{$split1[10]}=$temp;
    }
   }
  }
 }
}
close file1;

open file1, "<$gencc";
my (%gencc,%genccid);
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if (exists $split1[8]){
  unless ($split1[8] eq "Refuted Evidence" || $split1[8] eq "No Known Disease Relationship"){
   my $moi=$split1[10];
   if (exists $moi{$split1[10]}){
    $moi=$moi{$split1[10]};
   }
   if ($moi=~/AD/ || $moi=~/XLD/ || $moi=~/dominant/){
    $domgene{$split1[2]}=1;
   }

   $genccid{$split1[2]}=$split1[1];
   if (exists $gencc{$split1[2]}){
    $gencc{$split1[2]}="$gencc{$split1[2]}<br/>$split1[4] ($moi), $split1[8]";
   }
   else{
    $gencc{$split1[2]}="$split1[4] ($moi), $split1[8]";
   }
  }
 }
}
close file1;

open file1, "<$goterms";
my %goterms;
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 my @split4=split /\|/,$split1[1];
 my $temp=$split4[0];
 for (my $i=1; $i<=$#split4; $i++){
  $temp="$temp, $split4[$i]";
 }
 $split1[1]=$temp;

 $goterms{$split1[0]}=$split1[1];
}
close file1;

open file1, "<$expression";
my %heartexp;
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 $heartexp{$split1[0]}=$split1[1];
}
close file1;

open file1, "<$coding_exons";
open out1, ">$dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script07_file1.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 if ($split1[1] eq ""){
  $split1[1]=$split1[0];
 }
 if ($split1[2] eq $chr && $split1[4]>=$start && $split1[3]<=$end){
  my $coding_in=0;
  my @split14=split /\s/,$split1[5];
  LOOP1: foreach my $split14 (@split14){
   my @split15=split /\|/,$split14;
   if ($split15[2]>=$start && $split15[1]<=$end){
    $coding_in=1; last LOOP1;
   }
  }

  my (@split7,@split8,@split9,@split10,@split11);
  if (exists $enst{$split1[0]}){
   @split7=split /\s/,$enst{$split1[0]};
   @split8=split /\s/,$exon_ct{$split1[0]};
   @split9=split /\s/,$coding_exon{$split1[0]};
   @split10=split /\s/,$all_exon{$split1[0]};
   @split11=split /\s/,$mane{$split1[0]};
  }

  if ($coding_in==1){
   my $overlap="Partial";
   if ($start<=$split1[3] && $end>=$split1[4]){
    $overlap="Full";
   }
   if ($split1[3]<=$start && $split1[4]>=$end){
    $overlap="Intragenic";
   }

   my $trans="";
   if ($#split7>=0){
    for (my $i=0; $i<=$#split7; $i++){
     my $part=""; my $fullfrom=""; my $fullto=""; my $affected_bp=0; my $tot_bp=0;
     my @split16=split /\|\|/,$split9[$i];
     my @split17=split /\|\|/,$split10[$i];
     for (my $j=0; $j<=$#split17; $j++){
      my $exonnum=$j+1;
      my @split18=split /\:/,$split17[$j];
      if ($start<=$split18[0] && $end>=$split18[1]){
       if ($fullfrom eq ""){
        $fullfrom=$exonnum; $fullto=$exonnum;
       }
       elsif ($exonnum<$fullfrom){$fullfrom=$exonnum}
       elsif ($exonnum>$fullto){$fullto=$exonnum}
      }
      elsif ($end>=$split18[0] && $start<=$split18[1]){
       if ($part eq ""){$part=$exonnum}
       else{$part="$part, $exonnum"}
      }
     }

     for (my $j=0; $j<=$#split16; $j++){
      my @split19=split /\:/,$split16[$j];
      $tot_bp+=($split19[1]-$split19[0]+1);
      if ($end>=$split19[0] && $start<=$split19[1]){
       my @split20=@split19; push @split20, $start; push @split20, $end;
       @split20=sort{$a<=>$b}@split20;
       $affected_bp+=($split20[2]-$split20[1]+1);
      }
     }

     my $full="";
     if ($fullfrom ne ""){
      if ($fullfrom eq $fullto){
       $full="# $fullfrom";
      }
      else{
       $full="# $fullfrom-$fullto";
      }
     }
     else{
      $full="-";
     }

     if ($part eq ""){
      $part="-";
     } 
     else{
      $part="# $part";
     }

     my $trans_out="<tr>";
     if ($split11[$i]==1){
      $trans_out="$trans_out<td><a href='https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=$split1[0];t=$split7[$i]' target='_blank'>$split7[$i]</a> <span class='label label-warning'>MANE</span></td>";
     }
     else{
      $trans_out="$trans_out<td><a href='https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=$split1[0];t=$split7[$i]' target='_blank'>$split7[$i]</a></td>";
     }
     $trans_out="$trans_out<td>$split8[$i]</td><td>$full</td><td>$part</td>";
     if ($affected_bp%3==0){
      my $rt_cds=int($affected_bp*1e3/$tot_bp)/10;
      $trans_out="$trans_out<td><span class='label label-info'>In-frame</span> $affected_bp bp ($rt_cds%)</td>";
     }
     else{
      $trans_out="$trans_out<td><span class='label label-danger'>Out-of-frame</span> $affected_bp bp </td>";
     }
     $trans_out="$trans_out</tr>";
     $trans="$trans$trans_out";
    }
   }
   if ($trans eq ""){
    $trans="Transcript information not available for this gene. Please check Ensembl, UCSC or gnomAD manually.";
   }
   else{
    $trans="<table border='1px' bordercolor='#bdbdbd'><tr><th>Transcript (Sorted by MANE & <a href='https://gnomad.broadinstitute.org/gene/$split1[0]?dataset=gnomad_r2_1' target='_blank'>Expression</a>)&emsp;&emsp;</th><th>Total Exon(s)&emsp;</th><th>Exon Number(s) - Full Overlap&emsp;</th><th>Exon Number(s) - Partial Overlap&emsp;</th><th>Coding Sequence Overlap&emsp;</th></tr>$trans";
   }
   $trans="$trans</table><hr><iframe src='https://gnomad.broadinstitute.org/gene/$split1[0]?dataset=gnomad_r2_1' height='500px' width='95%'></iframe>";

   my $pli="NA"; my $loeuf="NA"; my $prec="NA"; my $pnull="NA";
   my $dom=0; my $omim="NA"; my $gencc="NA";
   my $heartexp="NA"; my $goterms="NA"; my $hi="NA"; my $ts="NA";
   my $synoposis="NA";

   if (exists $pli{$split1[0]}){
    if ($pli{$split1[0]} ne "NA"){
     $pli=int(1e2*$pli{$split1[0]})/1e2;
    }
    if ($loeuf{$split1[0]} ne "NA"){
     $loeuf=int(1e2*$loeuf{$split1[0]})/1e2;
    }
    if ($prec{$split1[0]} ne "NA"){
     $prec=int(1e2*$prec{$split1[0]})/1e2;
    } 
    if ($pnull{$split1[0]} ne "NA"){
     $pnull=int(1e2*$pnull{$split1[0]})/1e2;
    }
   }

   if (exists $domgene{$split1[0]} || exists $domgene{$split1[1]}){
    $dom=1;
   }

   if (exists $goterms{$split1[1]}){
    $goterms=$goterms{$split1[1]};
   }
 
   if (exists $heartexp{$split1[1]}){
    $heartexp=$heartexp{$split1[1]};
   }

   if (exists $hi{$split1[1]}){
    $hi="<a href='https://dosage.clinicalgenome.org/clingen_gene.cgi?sym=$split1[1]' target='_blank'>$hi{$split1[1]}</a>";
   }
 
   if (exists $ts{$split1[1]}){
    $ts="<a href='https://dosage.clinicalgenome.org/clingen_gene.cgi?sym=$split1[1]' target='_blank'>$ts{$split1[1]}</a>";
   }
   if ($type eq "DEL"){
    $ts="NA";
   }

   if (exists $omim{$split1[0]}){
    my @split2=split /\::/,$omim{$split1[0]};
    my $omim_out="";
    foreach my $split2 (@split2){
     my @split3=split /\|/,$split2;
     if ($omim_out eq ""){
      $omim_out="<a href='https://www.omim.org/entry/$split3[1]' target='_blank'>$split3[0]</a>";
     }
     else{
      $omim_out="$omim_out<br/><a href='https://www.omim.org/entry/$split3[1]' target='_blank'>$split3[0]</a>";
     } 
 
     if ($synoposis eq "NA"){
      $synoposis=$split3[1];
     }
     else{
      $synoposis="$synoposis,$split3[1]";
     }
    }
    $omim=$omim_out;
   }
 
   if (exists $gencc{$split1[1]}){
    $gencc="<a href='https://search.thegencc.org/genes/$genccid{$split1[1]}' target='_blank'>$gencc{$split1[1]}</a>";
   }

   my $strand_out="Unknown";
   if (exists $strand{$split1[0]}){
    $strand_out=$strand{$split1[0]};
   }
 
   if ($dom==1){
    print out1 "<a href='https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=$split1[0]' target='_blank'>$split1[1]</a><br/><span class='label label-warning'>Dominant</span>";
   }
   else{
    print out1 "<a href='https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=$split1[0]' target='_blank'>$split1[1]</a>";
   }
   if (exists $noncoding_exc{$split1[1]}){
    print out1 "<br/><span class='label label-warning'>Noncoding Pathogenic</span>\t";
   } 
   else{
    print out1 "\t";
   }
   
   print out1 "$split1[2]\t$split1[3]\t$split1[4]\t$overlap\t$hi\t$ts\t$pli\t$loeuf\t$prec\t$pnull\t$heartexp\t$goterms\t$omim\t$gencc\t";
   if ($synoposis ne "NA"){
    print out1 "<a href='https://www.omim.org/clinicalSynopsis/table?mimNumber=$synoposis' target='_blank' class='label label-default'>Clinical Synoposis</a></br>";
    if (exists $omimgene{$split1[0]}){
     print out1 "<a href='https://www.omim.org/allelicVariants/$omimgene{$split1[0]}' target='_blank' class='label label-default'>Allelic Variants</a>";
    }
    print out1 "</br></br>";
   }
   print out1 "<a href='https://gtexportal.org/home/gene/$split1[1]' target='_blank' class='label label-default'>GTEx</a><br/>";
   print out1 "<a href='https://gnomad.broadinstitute.org/gene/$split1[0]?dataset=gnomad_r2_1' target='_blank' class='label label-default'>gnomAD-Gene</a><br/>";
   print out1 "<a href='https://www.ncbi.nlm.nih.gov/clinvar/?term=$split1[1]%5Bgene%5D' target='_blank' class='label label-default'>ClinVar</a><br/>";
   print out1 "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=$split1[1]' target='_blank' class='label label-default'>GeneCards</a><br/>";
   print out1 "\t";
   open file2, "<$search";
   while (<file2>){
    chomp;
    my @split5=split /\t/,$_;
    my @split6=split /\s/,$split5[0];
    print out1 "<a href='https://www.google.com/search?q=%22$split1[1]%22";
    foreach my $split6 (@split6){
     print out1 "+$split6";
    }
    print out1 "' target='_blank' class='label label-default'>$split5[1]</a>   ";
   }
   print out1 "\t$dom\t$strand_out\t$trans\n";
  }
 }
}
close file1;
close out1;
system ("gzip -f $dir/app_temp_file/$proband/$proband.$chr.$start.$end.$type.script07_file1.txt");

exit 2;
