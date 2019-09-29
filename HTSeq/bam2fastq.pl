# convert paired bam file to two fastq files
# prerequisite: samtools (>=1.4), sambamba
# bam_files: one bam file including path
# output: the output folder. Two fastq files will be written here: fastq1.fastq and fastq2.fastq. You can further gzip them
# thread: number of threads to use
#!/usr/bin/perl
use strict;
use warnings;

my ($bam_files,$output,$thread)=@ARGV;
my ($header,$line,%single,$bam_file);

open(FILE_OUT1,">".$output."/fastq1.fastq") or die "Error: cannot write to ".$output."/fastq1.fastq!\n";
open(FILE_OUT2,">".$output."/fastq2.fastq") or die "Error: cannot write to ".$output."/fastq2.fastq!\n";

foreach $bam_file (split(",",$bam_files))
{
  # bam2fq
  print "bam2fq for ".$bam_file."\n";
  system("sambamba sort -o ".$output."/sorted.bam -N -t ".$thread." ".$bam_file);
  system("samtools bam2fq -1 ".$output."/tmp1.fq -2 ".$output."/tmp2.fq -n --threads ".$thread." ".$output."/sorted.bam");
  unlink($output."/sorted.bam");

  # identify singletons
  %single=();
  open(FILE_IN1,$output."/tmp1.fq");
  open(FILE_IN2,$output."/tmp2.fq");

  while ($line=<FILE_IN1>) 
  {
    add(\%single,$line); 
    <FILE_IN1>;<FILE_IN1>;<FILE_IN1>;

    $line=<FILE_IN2>;
    if (defined $line)
    {
      add(\%single,$line);
      <FILE_IN2>;<FILE_IN2>;<FILE_IN2>;
    }
  }

  while ($line=<FILE_IN2>)
  {
    add(\%single,$line);
    <FILE_IN2>;<FILE_IN2>;<FILE_IN2>;
  }

  close(FILE_IN1);
  close(FILE_IN2);

  # print paired reads
  open(FILE_IN,$output."/tmp1.fq");
  while ($header=<FILE_IN>) 
  {
    $line=<FILE_IN>;
    $line.=<FILE_IN>;
    $line.=<FILE_IN>;
    if (exists $single{$header}) {delete $single{$header}} else {print FILE_OUT1 $header.$line;}
  }
  close(FILE_IN);

  open(FILE_IN,$output."/tmp2.fq");
  while ($header=<FILE_IN>)
  {
    $line=<FILE_IN>;
    $line.=<FILE_IN>;
    $line.=<FILE_IN>;
    if (exists $single{$header}) {delete $single{$header}} else {print FILE_OUT2 $header.$line;}
  }
  close(FILE_IN);

  unlink($output."/tmp1.fq");
  unlink($output."/tmp2.fq");
}

close(FILE_OUT1);
close(FILE_OUT2);

sub add
{
  my ($single_ref,$line)=@_;
  if (exists $single_ref->{$line})
  {
    delete($$single_ref{$line});
  }else
  {
    $$single_ref{$line}=1;
  }
}

