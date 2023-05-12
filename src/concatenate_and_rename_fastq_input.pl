#!/usr/bin/perl

my ($data_dir, $targets_file, $sampleID) = ($ARGV[0]."/", $ARGV[1], $ARGV[2]);


# get analysis ID from the targets file

open FPIN, "<".$targets_file or die;
my @IDS;

while (<FPIN>) {

	push @IDS, $1 if (/^([^\,]+)\,$sampleID[\,\n\r\f\t]/); 

}
close FPIN;

die "too many instances of $sampleID in targets file" if (scalar @IDS > 1);
die "no instances of $sampleID in targets file" if (scalar @IDS == 0);


my $id = @IDS[0];

opendir DIR, $data_dir or die;
my @files = grep { /^$id\_?.*\_R1_/ || /^$id\_1/ } readdir DIR;
closedir DIR;


my $string = "cat";
foreach my $file (@files) {

	$string .= " ".$data_dir.$file;
}

$string .= " > ".$sampleID."/".$sampleID."_1.fastq.gz\n";
`mkdir $sampleID`;
`$string`;

$string =~ s/_R1_/_R2_/g;
$string =~ s/_1.fastq.gz/_2.fastq.gz/g;

`$string`;



