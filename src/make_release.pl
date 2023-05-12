#!/usr/bin/perl

# ./make_release.pl /scratch/c.mtera/nextflow/nf003 nf003 targets.csv /scratch/c.mtera/data.nextflow/ds0140  output resources Mus_musculus.GRCm38.98.gtf
my ($work_dir, $study_id, $sample_sheet, $data_dir, $output_dir, $resources_dir, $gtf_file) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4], $ARGV[5], $ARGV[6]);

my @sample_array;

# get sample names

open FPIN, "<".$data_dir."/".$sample_sheet or die;
my $line = <FPIN>;

while (<FPIN>) {
	die unless (/^[^\,]+\,([^\,\f\r\n]+)[\,\f\r\n]+/);
	push @sample_array, $1;
}
close FPIN;

my $ds;

my $rel_dir = "rds-release-".$study_id."-".`date +%Y%m%d`;
$rel_dir =~ s/\n//;

# get fastq files

`mkdir -p ${rel_dir}/fastq`;


foreach my $sample (@sample_array) {


	if (-e $work_dir."/".$output_dir."/".$sample."/".$sample."_1.fastq.gz") {
		`cp ${work_dir}/${output_dir}/${sample}/${sample}_1.fastq.gz $rel_dir/fastq`;
	} elsif (-e $work_dir."/output/".$sample."/".$sample."_1.fastq.gz") {
		`cp ${work_dir}/output/${sample}/${sample}_1.fastq.gz $rel_dir/fastq`;
	} else {
		die "input: ".$sample."_1.fastq.gz does not exist";
	}

	if (-e $work_dir."/".$output_dir."/".$sample."/".$sample."_2.fastq.gz") {
                `cp ${work_dir}/${output_dir}/${sample}/${sample}_2.fastq.gz $rel_dir/fastq`;
        } elsif (-e $work_dir."/output/".$sample."/".$sample."_2.fastq.gz") {
                `cp ${work_dir}/output/${sample}/${sample}_2.fastq.gz $rel_dir/fastq`;
        } else {
                die "input: ".$sample."_2.fastq.gz does not exist";
        }

}

# get bams and bai and gtf for igv

#`mkdir -p ${rel_dir}/igv/`;

#foreach my $sample (@sample_array) {

#        if (-e $work_dir."/".$output_dir."/".$sample."/".$sample.".markdup.bam") {
#                `cp ${work_dir}/${output_dir}/${sample}/${sample}.markdup.bam $rel_dir/igv`;
#        } elsif (-e $work_dir."/output/".$sample."/".$sample.".markdup.bam") {
#                `cp ${work_dir}/output/${sample}/${sample}.markdup.bam $rel_dir/igv`;
#        } else {
#                die "input: ".$sample.".markdup.bam does not exist";
#        }

#        if (-e $work_dir."/".$output_dir."/".$sample."/".$sample.".markdup.bam.bai") {
#                `cp ${work_dir}/${output_dir}/${sample}/${sample}.markdup.bam.bai $rel_dir/igv`;
#        } elsif (-e $work_dir."/output/".$sample."/".$sample.".markdup.bam.bai") {
#                `cp ${work_dir}/output/${sample}/${sample}.markdup.bam.bai $rel_dir/igv`;
#        } else {
#                die "input: ".$sample." does not exist";
#        }

#        if (-e $work_dir."/".$output_dir."/".$sample."/".$sample.".rmdup.bam") {
#                `cp ${work_dir}/${output_dir}/${sample}/${sample}.rmdup.bam $rel_dir/igv`;
#        } elsif (-e $work_dir."/output/".$sample."/".$sample.".rmdup.bam") {
#                `cp ${work_dir}/output/${sample}/${sample}.rmdup.bam $rel_dir/igv`;
#        } else {
#                die "input: ".$sample." does not exist";
#        }

#        if (-e $work_dir."/".$output_dir."/".$sample."/".$sample.".rmdup.bam.bai") {
#                `cp ${work_dir}/${output_dir}/${sample}/${sample}.rmdup.bam.bai $rel_dir/igv`;
#        } elsif (-e $work_dir."/output/".$sample."/".$sample.".rmdup.bam.bai") {
#                `cp ${work_dir}/output/${sample}/${sample}.rmdup.bam.bai $rel_dir/igv`;
#        } else {
#                die "input: ".$sample." does not exist";
#        }
#}

#die "resources: ".$gtf_file." does not exist" unless (-e $work_dir."/".$resources_dir."/".$gtf_file);
#`cp ${work_dir}/${resources_dir}/${gtf_file} $rel_dir/igv`;

# get bamtools and multiqc reports 

`mkdir -p ${rel_dir}/qc/`;

die "bamtools: all_mapping_stats.txt does not exist" unless (-e $work_dir."/".$output_dir."/all_mapping_stats.txt");
`cp ${work_dir}/${output_dir}/all_mapping_stats.txt $rel_dir/qc`;

die "multiqc: multiQC.html does not exist" unless (-e $work_dir."/".$output_dir."/multiQC.html");
`cp ${work_dir}/${output_dir}/multiQC.html $rel_dir/qc`;


# make resources

`mkdir -p ${rel_dir}/resources/`;

`cp ${work_dir}/nextflow.config $rel_dir/resources/`;
`cp ${work_dir}/main.nf $rel_dir/resources/`;
`cp ${data_dir}/${sample_sheet} $rel_dir/resources/`;
`cp ${work_dir}/${output_dir}/all*txt $rel_dir/resources/`;
`cp ${work_dir}/${output_dir}/multiQC.html $rel_dir/resources/`;
`cp -r ${work_dir}/src $rel_dir/resources/`;
`cp ${work_dir}/slurm/* $rel_dir/resources/`;
`cp ~/bin/run_nextflow.sh $rel_dir/resources/`;

# make our release

my $rel_dir = "onedrive-release-".$study_id."-".`date +%Y%m%d`;
$rel_dir =~ s/\n//;
`mkdir -p ${rel_dir}`;

`cp ${work_dir}/nextflow.config $rel_dir`;
`cp ${work_dir}/main.nf $rel_dir`;
`cp ${data_dir}/${sample_sheet} $rel_dir`;
`cp ${work_dir}/${output_dir}/all*txt $rel_dir`;
`cp ${work_dir}/${output_dir}/multiQC.html $rel_dir`;
`cp -r ${work_dir}/src $rel_dir`;
`cp ${work_dir}/slurm/* $rel_dir`;
`cp ~/bin/run_nextflow.sh $rel_dir`;







