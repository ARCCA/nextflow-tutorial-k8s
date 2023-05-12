#!/usr/bin/perl

#suppliedID,analysisID
#R122-H-001,wt1
#R122-H-002,wt2
#R122-H-003,wt3

my ($sample_sheet, $work_dir, $output_dir, $link_dir) = ($ARGV[0], $ARGV[1]."/", $ARGV[2]."/", $ARGV[3]."/");


my @sample_array;

open FPIN, "<".$sample_sheet or die;
my $line = <FPIN>;

while (<FPIN>) {
	die unless (/^[^\,]+\,([^\,\f\r\n]+)[\,\f\r\n]/);
	push @sample_array, $1;
}
close FPIN;

my $ds;

foreach my $sample (@sample_array) {

print $sample."\n";
	# total starting reads pairs
	
     	open FPIN, "<".$link_dir."/".$sample."/".$sample."_1_fastqc.html" or die "[".$work_dir."][".$link_dir."]/[".$sample."]/[".$sample."]_1_fastqc.html";

        while (<FPIN>) {

                if (/<td>Total Sequences<\/td><td>(\d+)<\/td>/) {
                        $ds->{$sample}->{"total"} = $1;
                        $ds->{$sample}->{"ptotal"} = "100";

                }
        }

        close FPIN;

        # total read pairs after trimming
       
	open FPIN, "<".$link_dir."/".$sample."/".$sample.".trimmed_1_fastqc.html" or die $work_dir."/".$output_dir."/".$sample."/".$sample.".trimmed_1_fastqc.html"; 
        
        while (<FPIN>) {
        
       		if (/<td>Total Sequences<\/td><td>(\d+)<\/td>/) {
                	$ds->{$sample}->{"trim"} = $1;
                        $ds->{$sample}->{"ptrim"} = (($1 / $ds->{$sample}->{"total"}) * 100);
		}
	}
      
       close FPIN;


      open FPIN, "<".$link_dir."/".$sample."/".$sample.".markdup.stats.txt" or die $work_dir."/".$output_dir."/".$sample."/".$sample.".markdup.stats.txt";

        while (<FPIN>) {

                if (/^Mapped reads:\s+(\d+)[^\(]+\(([^\%]+)[\%]/) {
                        $ds->{$sample}->{"mapped"} = ($1 / 2);
                        $ds->{$sample}->{"pmapped"} = (($ds->{$sample}->{"mapped"} / $ds->{$sample}->{"total"} * 100));
                } elsif (/^Forward strand:\s+(\d+)[^\(]+\(([^\%]+)[\%]/) {
                        $ds->{$sample}->{"forward"} = ($1 / 2);
                        $ds->{$sample}->{"pforward"} = $2;
                } elsif (/^Reverse strand:\s+(\d+)[^\(]+\(([^\%]+)[\%]/) {
                        $ds->{$sample}->{"reverse"} = ($1 / 2);
                        $ds->{$sample}->{"preverse"} = $2;
                } elsif (/^Duplicates:\s+(\d+)[^\(]+\(([^\%]+)[\%]/) {
                        $ds->{$sample}->{"duplicate"} = ($1 / 2);
	             $ds->{$sample}->{"pduplicate"} = $2;
               }
        }

    close FPIN;

  # open picard file

    open FPIN, "<".$link_dir."/".$sample."/".$sample.".markdup.picardstats.txt" or die $work_dir."/".$output_dir."/".$sample."/".$sample.".markdup.picardstats.txt";

        while (<FPIN>) {

                if (/^MEDIAN_INSERT_SIZE/) {
			my $line = <FPIN>;
			die unless ($line =~ /^([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t/);

			$ds->{$sample}->{"median_insert_size"} = $1;
			$ds->{$sample}->{"mode_insert_size"} = $2;
			$ds->{$sample}->{"mean_insert_size"} = $3;
			$ds->{$sample}->{"standard_deviation"} = $4;
			last;
                } 
        }

    close FPIN;

}

open FPOUT, ">all_mapping_stats.txt" or die; 

print FPOUT "SampleID\t";
print FPOUT "No. read-pairs\t";
print FPOUT "No. read-pairs (\%)\t";
print FPOUT "No. read-pairs after trimming\t";
print FPOUT "No. read-pairs after trimming (\%)\t";
print FPOUT "No. read-pairs after mapping\t";
print FPOUT "No. read-pairs after mapping (\%)\t";
print FPOUT "No. read-pairs after mapping on forward strand\t";
print FPOUT "No. read-pairs after mapping on forward strand (\%)\t";
print FPOUT "No. read-pairs after mapping on reverse strand\t";
print FPOUT "No. read-pairs after mapping on reverse strand (\%)\t";
print FPOUT "No. read-pairs duplicated\t";
print FPOUT "No. read-pairs duplicated (\%)\t";
print FPOUT "Median insert length\t";
print FPOUT "Mode insert length\t";
print FPOUT "Mean insert length\t";
print FPOUT "Standard deviation of insert length\t";
print FPOUT "\n";

foreach my $sample (@sample_array) {


	print FPOUT $sample."\t";
	print FPOUT $ds->{$sample}->{"total"}."\t".$ds->{$sample}->{"ptotal"}."\t";
	print FPOUT $ds->{$sample}->{"trim"}."\t".$ds->{$sample}->{"ptrim"}."\t";
	print FPOUT $ds->{$sample}->{"mapped"}."\t".$ds->{$sample}->{"pmapped"}."\t";
	print FPOUT $ds->{$sample}->{"forward"}."\t".$ds->{$sample}->{"pforward"}."\t";
	print FPOUT $ds->{$sample}->{"reverse"}."\t".$ds->{$sample}->{"preverse"}."\t";
	print FPOUT $ds->{$sample}->{"duplicate"}."\t".$ds->{$sample}->{"pduplicate"}."\t";
	print FPOUT $ds->{$sample}->{"median_insert_size"}."\t".$ds->{$sample}->{"mode_insert_size"}."\t".$ds->{$sample}->{"mean_insert_size"}."\t".$ds->{$sample}->{"standard_deviation"}."\n";
}

close FPOUT;

