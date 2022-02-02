#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::Util qw(shuffle);
my %opts;
getopt ('n:',\%opts);
my $name=$opts{"n"};

my @arr1;
open file1, "<./user_data/$name.filtered.txt";
while (<file1>){
 chomp;
 push @arr1, $_;
}
close file1;

@arr1=shuffle(@arr1);
open out1, ">./user_data/$name.hg19.variant_list.txt";
foreach my $arr1 (@arr1){ 
 my @split1=split /\t/,$arr1;
 print out1 "$split1[6]\t$split1[0]\t$split1[1]\t$split1[2]\t$split1[3]\n";
}
close file1;
close out1;

exit 2;
