#!/usr/bin/perl
use strict;
use warnings;

#fastq1: first mate
#fastq2: second mate
#fc: featurecount files
#output_fastq: output fastq file
#n: number of reads to extract each time
my ($fastq1,$fastq2,$fc,$output_fastq,$n)=@ARGV;
my (@items,$name,$i,$j,%reads,$line11,$line12,$line13,$line14,$line21,$line22,$line23,$line24);

$i=0;
$name="";

open(FILE_FC,$fc);
open(FILE_FQ1,$fastq1);
open(FILE_FQ2,$fastq2);
open(FILE_UNMAP1,">".$output_fastq."_1.fq");
open(FILE_UNMAP2,">".$output_fastq."_2.fq");

while (<FILE_FC>)
{
  @items=split("\t",$_);
  if ($items[0] eq $name) # same RNA-Seq read name as above
  {
    if ($items[1] ne "Unassigned_NoFeatures") {$reads{$name}->{"status"}=0;} # assume Unassigned_NoFeatures is the only one we want. We can discuss this later
  }else # found a different one
  {
    if (exists $reads{$name} && $reads{$name}->{"fastq1_seq"} ne "") # process the previous read on featurecounts' side
    {
      if ($reads{$name}->{"status"}==1) # if matched and wanted 
      {
        print FILE_UNMAP1 new_fq($name,$reads{$name}->{"fastq1_seq"},$reads{$name}->{"fastq1_qual"});
        print FILE_UNMAP2 new_fq($name,$reads{$name}->{"fastq2_seq"},$reads{$name}->{"fastq2_qual"});
      }
      delete $reads{$name};
    }

    if ($i==$n) # process reads on fastq files' side
    {
      $i=0;
      foreach $j (0..$n-1) # read n reads one by one
      {
        $line11=<FILE_FQ1>;
        $line12=<FILE_FQ1>;
        $line13=<FILE_FQ1>;
        $line14=<FILE_FQ1>;

        $line21=<FILE_FQ2>;
        $line22=<FILE_FQ2>;
        $line23=<FILE_FQ2>;
        $line24=<FILE_FQ2>;

        @items=split(" ",$line11); # find read name
        $items[0]=~s/^@//;
        $name=$items[0];

        if (exists $reads{$name}) # we have seen the read before
        {
          if ($reads{$name}->{"status"}==1)
          {
            print FILE_UNMAP1 new_fq($name,$line12,$line14);
            print FILE_UNMAP2 new_fq($name,$line22,$line24);
          }
          delete $reads{$name};
        }else # if not, add it to the hash
        {
          $reads{$name}={status=>1,fastq1_seq=>$line12,fastq1_qual=>$line14,fastq2_seq=>$line22,fastq2_qual=>$line24};
        }
      }
    }

    # back to featurecounts' side
    $i++;
    $name=$items[0];
    if (! exists $reads{$name}) {$reads{$name}={status=>1,fastq1_seq=>"",fastq1_qual=>"",fastq2_seq=>"",fastq2_qual=>""};} # assume we want this read by default
    if ($items[1] ne "Unassigned_NoFeatures") {$reads{$name}->{"status"}=0;}
  } 
}

# process the last few reads on fastq files' side
while ($line11=<FILE_FQ1>) # read n reads one by one
{
  $line12=<FILE_FQ1>;
  $line13=<FILE_FQ1>;
  $line14=<FILE_FQ1>;

  $line21=<FILE_FQ2>;
  $line22=<FILE_FQ2>;
  $line23=<FILE_FQ2>;
  $line24=<FILE_FQ2>;

  @items=split(" ",$line11); # find read name
  $items[0]=~s/^@//;
  $name=$items[0];

  if (exists $reads{$name}) # we have seen the read before
  {
    if ($reads{$name}->{"status"}==1)
    {
      print FILE_UNMAP1 new_fq($name,$line12,$line14);
      print FILE_UNMAP2 new_fq($name,$line22,$line24);
    }
    delete $reads{$name};
  }else # if not, add it to the hash
  {
    $reads{$name}={status=>1,fastq1_seq=>$line12,fastq1_qual=>$line14,fastq2_seq=>$line22,fastq2_qual=>$line24};
  }
}
if ((scalar keys %reads)!=0) {print "Leftover reads in the hash! Something is not right!\n";}

close(FILE_UNMAP1);
close(FILE_UNMAP2);
close(FILE_FQ1);
close(FILE_FQ2);
close(FILE_FC);

sub new_fq
{
  my ($name,$seq,$qual)=@_;
  "@".$name."\n".$seq."+\n".$qual;
}


