#/usr/bin/perl
use strict;

#Define hash of SNPs
my %snps_hash;

#File containing SNPs to be added to the LDscore annotations
my $snps=$ARGV[0];

#Name of the annotation
my $name=$ARGV[1];

#Reads file with SNPs to add to annotation and store them in a hash
open my $snps_file, '<', $snps or die "What have you done?! $!";
while (<$snps_file>) {
  chomp;                
  $snps_hash{$_} = 1;  
}
close $snps_file;

#Open Annotation file of LDSC
my $ld_scores_annot=$ARGV[2];
open my $ld_scores_file, '<', $ld_scores_annot or die "Oh no, not again.. $!";

#Print header 
my $header = <$ld_scores_file>;
chomp $header;
print "CHR\tBP\tSNP\tCM\t". $name ."\n";

#if SNP in annotation is also in hash -> Add 1 in the new annotation, else add 0.
while (<$ld_scores_file>) {
  chomp;
  my ($chr,$bp, $snp,$cm,$base) = split '\t';
  my $line = $_;
  #print $line . "\n";
  if (exists($snps_hash{$snp})){
  	print "$chr\t$bp\t$snp\t$cm\t1\n";
  }
  else {
  	print "$chr\t$bp\t$snp\t$cm\t0\n";
  }	
}
close $ld_scores_file;