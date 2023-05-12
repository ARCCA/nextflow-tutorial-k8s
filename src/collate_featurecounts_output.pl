#!/usr/bin/perl

my ($sample_sheet, $work_dir, $featurecounts_dir, $resources_dir, $gtf_file, $link_dir) = ($ARGV[0], $ARGV[1]."/", $ARGV[2]."/", $ARGV[3]."/", $ARGV[4], $ARGV[5]);

my @sample_array;

open FPIN, "<".$sample_sheet or die;
my $line = <FPIN>;

while (<FPIN>) {
	die unless (/^[^\,]+\,([^\,]+)[\,\f\r\n]/);
	push @sample_array, $1;
}
close FPIN;

my $ds;

# get annotation from the GTF files

my $g;
my $t;

 open FPIN, "<".$work_dir."/".$resources_dir."/".$gtf_file or die $work_dir."/".$resources_dir."/".$gtf_file;
      while (<FPIN>) {
		if (/gene_id \"([^\"]+)\".+gene_name \"([^\"]+)\".+gene_biotype \"([^\"]+)\";/) {
			$g->{$1}->{"name"} = $2;
			$g->{$1}->{"type"} = $3;

		}

               if (/gene_id \"([^\"]+)\".+transcript_id \"([^\"]+)\".+gene_name \"([^\"]+)\".+gene_biotype \"([^\"]+)\";/) {
                        $t->{$2}->{"gene"} = $1;
                        $t->{$2}->{"name"} = $3;
                        $t->{$2}->{"type"} = $4;
                }
	}

close FPIN;


foreach my $feature ("genecount") {


	foreach my $mapping ("markdup", "rmdup") {

		my $ids;
		my $ds;
		my @gene_order;

		open FPIN, "<".$link_dir."/".$sample_array[0]."/".$sample_array[0].".".$mapping.".".$feature.".txt" or die "<".$work_dir."/".$link_dir."/".$sample_array[0]."/".$sample_array[0].".".$mapping.".".$feature.".txt";
		<FPIN>;
		<FPIN>;
		while (<FPIN>) {

			die $_ unless (/^([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t\n]+)\n$/);

			$ids->{$1} = $2;
			push @gene_order, $1;

		}
		close FPIN;


		foreach my $sample (@sample_array) {

	               open FPIN, "<".$link_dir."/".$sample."/".$sample.".".$mapping.".".$feature.".txt" or die;
        	        <FPIN>;
                	<FPIN>;
                	while (<FPIN>) {

                        	die unless (/^([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t\n]+)\n$/);

                        	$ds->{$1}->{$sample} = $3;

                	}
                	close FPIN;

		}	

		open FPOUT, ">all.".$mapping.".".$feature.".txt" or die; 

		print FPOUT "ensemblGeneID\tgeneLength\tgeneName\tgeneBiotype";		

		foreach my $sample (@sample_array) {
			print FPOUT "\t".$sample;
		}

		print FPOUT "\n";

		foreach my $gene (@gene_order) {

			print FPOUT $gene;
			print FPOUT "\t".$ids->{$gene};
			print FPOUT "\t".$g->{$gene}->{"name"};
			print FPOUT "\t".$g->{$gene}->{"type"};
	
			foreach my $sample (@sample_array) {
				print FPOUT "\t".$ds->{$gene}->{$sample};
			}

			print FPOUT "\n";
		}

		close FPOUT;


	}
}
	
foreach my $feature ("transcriptcount") {

	foreach my $mapping ("markdup", "rmdup") {

		my $ids;
		my $ds;
		my @gene_order;

		open FPIN, "<".$link_dir."/".$sample_array[0]."/".$sample_array[0].".".$mapping.".".$feature.".txt" or die;
		<FPIN>;
		<FPIN>;
		while (<FPIN>) {

			die $_ unless (/^([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t\n]+)\n$/);

			$ids->{$1} = $2;
			push @gene_order, $1;

		}
		close FPIN;


		foreach my $sample (@sample_array) {

	               open FPIN, "<".$link_dir."/".$sample."/".$sample.".".$mapping.".".$feature.".txt" or die;
        	        <FPIN>;
                	<FPIN>;
                	while (<FPIN>) {

                        	die unless (/^([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t\n]+)\n$/);

                        	$ds->{$1}->{$sample} = $3;

                	}
                	close FPIN;

		}

		open FPOUT, ">all.".$mapping.".".$feature.".txt" or die;

		print FPOUT "ensemblTranscriptID\ttranscriptLength\tensemblGeneID\tgeneName\tgeneBiotype";

		foreach my $sample (@sample_array) {
			print FPOUT "\t".$sample;
		}

		print FPOUT "\n";

		foreach my $gene (@gene_order) {

			print FPOUT $gene;
			print FPOUT "\t".$ids->{$gene};
			print FPOUT "\t".$t->{$gene}->{"gene"};
			print FPOUT "\t".$t->{$gene}->{"name"};
			print FPOUT "\t".$t->{$gene}->{"type"};

			foreach my $sample (@sample_array) {
				print FPOUT "\t".$ds->{$gene}->{$sample};
			}

			print FPOUT "\n";
		}

		close FPOUT;


	}
}	
