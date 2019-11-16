#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long qw(:config posix_default no_ignore_case);
use Pod::Usage;
#use File::Tee qw(tee);
use File::Basename;
use Config;
use FindBin qw($Bin);
use threads;


my $target_scaf_file = '';     # -t
my $query_contig_file = '';    # -q
my $align_file = '';            # -a
my $align_file_Q = '';          # -aq
my $query_index = '';          # -qi
my @read_files = ();        # -r
my $read_format = 'fastq';  # -f
my $base_call_type = 'auto';    # -bq
my $read_length = 0;        # -l
my $insert_size = 400;      # -i
my $SD_insert = 40;         # -d
my @sam_target_align_file = (); # -st
my @sam_query_align_file = ();  # -sq
my $sam_dir = '';           # -sd
my @blast_target_align_file = ();
my @blast_query_align_file = ();
my $out_prefix = '';        # -p
my $min_match_length = 300;  # -mm
my $min_match_len_local = 20; # -ml
my $min_identity = 99;      # -mi
my $min_gap_size = 1;       # -mg
my $max_indel_size = 70;    # -is
my $min_subcon_size = 150;  #
my $max_query_singletone_cov = 60;  # -qsc
my $extend = 0;             # -et
my $LR = 0;                 # -lr
my $LR_coverage = 'auto';   # -lc
my $min_qalign_num = 1;     # -mq
my $connect_subcon = 0;     # -c
my $blast = 0;              # -b
my $adjust_score = 0;
my $heterozygosity = 0;     # -ht
my $thread_num = 5;         # -n
my $thread_connect = 0;     # -nc
my $iteration = 3;          # -it
my $ite_num = 1;            # -in
my $help;                   # -h

my $PE_align_tag = 0;
my @options = @ARGV;

GetOptions(
    'target_scaf|t=s' => \$target_scaf_file,
    'query_seq|q=s' => \$query_contig_file,
    'align_file|a=s' => \$align_file,
    'read_file|r=s{,}' => \@read_files,
    'query_index|qi=s' => \$query_index,
    'read_format|f' => \$read_format,
    'base_qual|bq=s' => \$base_call_type,
    'read_len|l=i' => \$read_length,
    'insert|i=i' => \$insert_size,
    'sd_insert|d=i' => \$SD_insert,
    'sam_talign|st=s{,}' => \@sam_target_align_file,
    'sam_qalign|sq=s{,}' => \@sam_query_align_file,
    'sam_dir|sd=s' => \$sam_dir,
    'prefix_out|p=s' => \$out_prefix,
    'min_match_len|mm=i' => \$min_match_length,
    'min_len_local|ml=i' => \$min_match_len_local,
    'min_identity|mi=i' => \$min_identity,
    'min_subcon|ms=i' => \$min_subcon_size,
    'min_gap_size|mg=i' => \$min_gap_size,
    'max_indel|is=i' => \$max_indel_size,
    'max_qscov|qsc=i' => \$max_query_singletone_cov,
    'conect_subcon|c' => \$connect_subcon,
    'blast|b' => \$blast,
    'extend|et' => \$extend,
    'ad_score|as=f' => \$adjust_score,
    'long_read|LR|lr' => \$LR,
    'lr_cov|lc=i' => \$LR_coverage,
    'min_qalign|mq=i' => \$min_qalign_num,
    'alignq|aq=s' => \$align_file_Q,
    'iterate|it=i' => \$iteration,
    'ite_num|in=i' => \$ite_num,
    'hetero|ht' => \$heterozygosity,
    'thread|n=i' => \$thread_num,
    'thread_connect|nc=i' => \$thread_connect,
    'help|h' => \$help
) or pod2usage(-verbose => 0);

pod2usage(-verbose => 0) if $help;
die "--target_scaf option is not specied: $!\n" if ($target_scaf_file eq '');
die "--query_contig option is not specied: $!\n" if ($query_contig_file eq '');
die "--prefix_out option is not specied: $!\n" if ($out_prefix eq '');
die "The string specified with --prefix_out option must not contain directory names but only a file prefix name: $!\n" if ($out_prefix =~ /\//);
die "--read_len option is not specied: $!\n" if ($read_length <= 0) and (@read_files > 0);
die "short read files are required when conducting with multiple numbers of iterations: $!\n" if ($iteration > 1) and (@read_files == 0) and (@sam_target_align_file > 0 or @sam_query_align_file > 0);
die "The number of read files should be even number (separate files of paired-end reads): $!\n" if (@read_files > 0) and (@read_files % 2 == 1);
die "The value specified with --insert should be read_len < insert < 20000: $!\n" if ($insert_size <= $read_length) or ($insert_size > 20000);

$min_identity = 95 if ($min_identity < 95);
$min_identity = 100 if ($min_identity > 100);
$min_match_length = 50 if ($min_match_length < 50);
$min_match_len_local = 15 if ($min_match_len_local < 15);
$min_gap_size = 1 if ($min_gap_size < 1);
if ($thread_connect == 0){
    $thread_connect = $thread_num;
}

if (($read_length >= 90) and ($insert_size < 350)){
    $insert_size -= $read_length - 75;
    $read_length = 75;
}

$extend = 0;

my $ite_dir = '';
if ($iteration > 1){
    $ite_dir = "iteration-$ite_num";
    system ("mkdir $ite_dir") unless (-d $ite_dir);
    $ite_dir = "iteration-$ite_num/";
    if ($ite_num > 1){
        my $pre_ite_num = $ite_num - 1;
        $target_scaf_file = "iteration-$pre_ite_num/$out_prefix.gapclosed.fa";
        $align_file = '';
        $align_file_Q = '';
        @sam_target_align_file = ();
        @sam_query_align_file = ();
    }
}

my $temp_dir = 'temp';
system ("mkdir $temp_dir") unless (-d $temp_dir);

my $blast_dir = 'db';

my $target_file = $target_scaf_file;
my $query_file = $query_contig_file;

my $prefix_target;
my $prefix_query;
if ($target_file =~ /-merged\./){
    $prefix_target = $1 if ($target_file =~ /([^\/]+)-merged\.\w+$/);
    $prefix_target = $1 if ($target_file =~ /([^\/]+)-merged$/) and ($target_file !~ /([^\/]+)-merged\.\w+$/);
}
else{
    $prefix_target = $1 if ($target_file =~ /([^\/]+)\.\w+$/);
    $prefix_target = $1 if ($target_file =~ /([^\/]+)$/) and ($target_file !~ /([^\/]+)\.\w+$/);
}
if ($query_file =~ /-merged\./){
    $prefix_query = $1 if ($query_file =~ /([^\/]+)-merged\.\w+$/);
    $prefix_query = $1 if ($query_file =~ /([^\/]+)-merged$/) and ($query_file !~ /([^\/]+)-merged\.\w+$/);
}
else{
    $prefix_query = $1 if ($query_file =~ /([^\/]+)\.\w+$/);
    $prefix_query = $1 if ($query_file =~ /([^\/]+)$/) and ($query_file !~ /([^\/]+)\.\w+$/);
}


#open(my $log_fh, '>', "$ite_dir$out_prefix.log");
#tee(STDERR, $log_fh);

my $log_file = "$out_prefix.submit.log";
open (OUTLOG, "> $log_file");
print OUTLOG "##### Submitted command #####\n";
pop @options;
pop @options;
print OUTLOG "gmcloser @options\n\n";

if ($ite_num == 1){
    print OUTLOG "##### Specified options #####\n [target file]: $target_file [query file]: $query_file\n [read file]: @read_files\n";
    print OUTLOG " [blast alignment file]: $align_file\n" if ($align_file ne '');
    print OUTLOG " [short read alignment file for target]: @sam_target_align_file\n" if (@sam_target_align_file > 0);
    print OUTLOG " [short read alignment file for query]: @sam_query_align_file\n" if (@sam_query_align_file > 0);
    print OUTLOG " [short read alignment file directory]: $sam_dir\n" if ($sam_dir ne '');
    print OUTLOG " [bowtie2 index file basename for query long read]: $query_index\n" if ($query_index ne '');
    print OUTLOG " [read length]: $read_length [insert]: $insert_size [sd_insert]: $SD_insert [min_gap]: $min_gap_size [max_indel]: $max_indel_size [hetero]: $heterozygosity [thread]: $thread_num\n";
    print OUTLOG " [min_match_len]: $min_match_length [min_identity]: $min_identity\n" if (@sam_target_align_file == 0) and (@sam_query_align_file == 0) and (@read_files == 0);
    print OUTLOG " [min_match_len_local]: $min_match_len_local [min_subcon_length] $min_subcon_size [max_query_singletone_cov]: $max_query_singletone_cov [adjust_score]: $adjust_score [connect_subcon]: $connect_subcon [extend_scaffold_termini]: $extend\n";
    print OUTLOG " [long_read]: $LR [long_read_coverage]: $LR_coverage [min_qalign]: $min_qalign_num [iteration]: $iteration\n" if ($LR == 1) and ($align_file_Q eq '');
    print OUTLOG " [long_read]: $LR [long_read_coverage]: $LR_coverage [min_qalign]: $min_qalign_num [iteration]: $iteration [query-query blast alignment file]: $align_file_Q\n" if ($LR == 1) and ($align_file_Q ne '');
    print OUTLOG " [prefix_out]: $out_prefix\n";
    
    print STDERR "##### Specified options #####\n [target file]: $target_file [query file]: $query_file\n [read file]: @read_files\n";
    print STDERR " [blast alignment file]: $align_file\n" if ($align_file ne '');
    print STDERR " [short read alignment file for target]: @sam_target_align_file\n" if (@sam_target_align_file > 0);
    print STDERR " [short read alignment file for query]: @sam_query_align_file\n" if (@sam_query_align_file > 0);
    print STDERR " [short read alignment file directory]: $sam_dir\n" if ($sam_dir ne '');
    print STDERR " [bowtie2 index file basename for query long read]: $query_index\n" if ($query_index ne '');
    print STDERR " [read length]: $read_length [insert]: $insert_size [sd_insert]: $SD_insert [min_gap]: $min_gap_size [max_indel]: $max_indel_size [hetero]: $heterozygosity [thread]: $thread_num\n";
    print STDERR " [min_match_len]: $min_match_length [min_identity]: $min_identity\n" if (@sam_target_align_file == 0) and (@sam_query_align_file == 0) and (@read_files == 0);
    print STDERR " [min_match_len_local]: $min_match_len_local [min_subcon_length] $min_subcon_size [max_query_singletone_cov]: $max_query_singletone_cov [adjust_score]: $adjust_score [connect_subcon]: $connect_subcon [extend_scaffold_termini]: $extend\n";
    print STDERR " [long_read]: $LR [long_read_coverage]: $LR_coverage [min_qalign]: $min_qalign_num [iteration]: $iteration\n" if ($LR == 1) and ($align_file_Q eq '');
    print STDERR " [long_read]: $LR [long_read_coverage]: $LR_coverage [min_qalign]: $min_qalign_num [iteration]: $iteration [query-query blast alignment file]: $align_file_Q\n" if ($LR == 1) and ($align_file_Q ne '');
    print STDERR " [prefix_out]: $out_prefix\n";
}

my $localtime = localtime;
print STDERR "\nJob start: $localtime\n" if ($iteration == 1);
print STDERR "\n##### Job start at iteration-$ite_num ##### : $localtime\n\n" if ($iteration > 1);
print OUTLOG "\nJob start: $localtime\n" if ($iteration == 1);
print OUTLOG "\n##### Job start at iteration-$ite_num ##### : $localtime\n\n" if ($iteration > 1);


=head1 SYNOPSIS

   GMcloser ver. 1.5
   
   Options:
   --target_scaf or -t <STR>       input target scaffold fasta file (e.g., scaf.fa) [mandatory]
   --query_seq or -q <STR>         input query contig (or long-read) fasta file (e.g., contig.fa) (if long reads are used, -lr option must be specified) [mandatory]
   --prefix_out or -p <STR>        prefix name of output files [mandatory]
   
   --read_file or -r <STR>         space-separated list of fastq or fasta files (or gzip compressed files) of paired-end reads (e.g., -r read_pe1.fq read_pe2.fq)
   --read_len or -l <INT>          read length of paired-end reads specified with the -r, -st, -sq, or -sd option (mean read length if read lengths are variable) [mandatory]
   --insert or -i <INT>            average insert size of paired-end reads [>read_len <20001, default: 400]
   --sd_insert or -d <INT>         standard deviation of insert size of paired-end reads [default: 40]
   --read_format or -f <STR>       fastq or fasta [default: fastq]
   --sam_talign or -st <STR>       space-separated list of sam alignment file(s) for target scaffolds, created in a single-end read alignment mode for paired-end reads (e.g., -sa tPE1.sam tPE2.sam, for paired-end read files PE1.fq and PE2.fq)
   --sam_qalign or -sq <STR>       space-separated list of sam alignment file(s) for query contigs, created in a single-end read alignment mode for paired-end reads (e.g., -sa qPE1.sam qPE2.sam, for paired-end read files PE1.fq and PE2.fq)
   --sam_dir or -sd <STR>          path of directory (i.e., bowtie_align) containing sam alignment files generated from a previous job with GMcloser (this can be used in place of -st and -sq option)
   --query_index or -qi <STR>      basename of bowtie2 index files for query long reads [optional]
   --align_file or -a <STR>        Nucmer alignment file for query against target [optional]

   --connect_subcon or -c          connect between gap-encompassing subcontig pairs with their original (not merged with query contigs) termini [default: false]
   --extend or -et                 extend scaffold temini with aligned query sequences [default: false] (When using long read query, this option is disabled in the current version)
   --blast or -b                   conduct alignment between target and query contigs with BLASTn [default: false] (Nucmer alignment by default)
   
   --min_match_len or -mm <INT>    minimum overlap length to be filtered for BLASTn alignments.
                                   Contig-alignments that satisfy both the values specified with -mm and -mi are selected, irrespective of any mapping rates of PE-reads. [>49, default: 300]
   --min_identity or -mi <INT>     minimum overlap identity (%) to be filtered for BLASTn alignments. Alignments are selected by combination with -mm option. [95~100, default: 99]
   --min_len_local or -ml <INT>    minimum overlap length for merging between neighbor subcontigs with YASS aligner [>14, default: 20]
   --min_subcon or -ms <INT>       minimum length of subcontigs, to be used for gapclosing [default: 150]
   --min_gap_size or -g <INT>      minimum length of gap, when spliting the target scaffold sequences into subcontigs [>0, default: 1]
   --max_indel or -is <INT>        maximum length of indel, observed in alignments between target subcontigs and query contigs. The alignments separated by the indel will be merged. [default: 70]
   --max_qsc or -qsc <INT>         maximum alignment coverage (%) of query singletones for target subcontigs (query with >= INT is excluded from query singletone output) [default: 60]
   --base_qual or -bq <STR>        base call quality format of fastq read file; illumina (phred64) or sanger (phred33) [default: auto]
   --ad_score or -as <FLOAT>       score to add to (subtract from) the standard threshold score for selection of correct contig-subcontig alignments (e.g., 1 or -1) [default: 0]
   --hetero or -ht                 heterozygosity factor (specify this if your sequenced genome is heterozygous (>0.2% difference of the haploid size)) [default: false]
   --thread or -n <INT>            number of threads (for machines with multiple processors), enabling all the alignment processes in parallel [default: 5]
   --thread_connect or -nc <INT>   number of threads (for machines with multiple processors), enabling the subcontig-connection process in parallel [default: number specified with --thread]
   --help or -h                    output help message
   
   [Options for long read query]
   --long_read or -LR              query sequence file is a fasta file of long reads (PacBio reads must be error-corrected) [default: false] (alignment was peformed with blast)
   --lr_cov or -lc <INT>           fold coverage of long reads for target scaffolds [default: auto ; automatically calculated by dividing a total length of query by a total length of target]
   --min_qalign or -mq <INT>       minimum number of queries that are aligned to either 5'- or 3'-terminus of a target subcontig [default: 1]
   --iterate or -it <INT>          number of iteration [default: 3]
   --alignq or -aq <STR>           BLASTn alignment file for query against query [optional]

=cut


my $target_name = '-';
my $query_name = '-';
my $Tseq = '';
my $Qseq = '';
my %target_seq;
my %target_Is;
my %target_Ie;
my %target_subcon_seq;
my %target_gap_size;
my %target_noscaf;
my %Tsubcon_info;
my $target_size = 0;

if (($sam_dir ne '') and (-d $sam_dir)){
    my @target_file_list = <$sam_dir/*-target-align-se-*>;
    my @query_file_list = <$sam_dir/*-query-align-se-*>;
    push @sam_target_align_file, @target_file_list;
    push @sam_query_align_file, @query_file_list;
}

$PE_align_tag = 1 if ((@read_files > 0) or ((@sam_target_align_file > 0) and (@sam_query_align_file > 0)));

# get target sequences

# get subcontig information (start and end positions for each scaffold) from a target scaffold fasta file

$target_name = '-';
my $total_target_subcon_len = 0;
open (FILE, $target_file) or die "$target_file is not found: $!\n";
while (my $line = <FILE>){
    $line =~ s/(\r\n|\n|\r)//g;
    if ($line =~ /^>/){
        if ($target_name ne '-'){
            $Tseq =~ s/^[Nn]+// if ($Tseq =~ /^[Nn]/);
            $Tseq =~ s/[Nn]+$// if ($Tseq =~ /[Nn]$/);
            $target_seq{$target_name} = uc $Tseq;
            $target_size += length $Tseq;
            if ($Tseq =~ /[Nn]{$min_gap_size}/){
                &assign_Tsubcont ($target_name, $Tseq);
            }
            else{
                my $count_subcon_04d = sprintf ("%04d", 1);
                my $subcon_name = $target_name . '/' . $count_subcon_04d;
                $target_subcon_seq{$subcon_name} = $Tseq;
                $target_gap_size{$subcon_name} = 0;
                $target_noscaf{$target_name} = 1;
                $total_target_subcon_len += length $Tseq;
                $target_Is{$target_name} = 1;
                $target_Ie{$target_name} = length $Tseq;
            }
            $Tseq = '';
        }
        $target_name = $1 if ($line =~ /^>(\S+)/);
        die "Taerget scaffold name must not contain '||' '==' or '/INT$', where 'INT' indicates integer and '$' the end of the string." if (($target_name =~ /\|\|/) or ($target_name =~ /==/) or ($target_name =~ /\/\d+$/));
    }
    else{
        $Tseq .= $line;
    }
}
$Tseq =~ s/^[Nn]+// if ($Tseq =~ /^[Nn]/);
$Tseq =~ s/[Nn]+$// if ($Tseq =~ /[Nn]$/);
$target_seq{$target_name} = uc $Tseq;
$target_size += length $Tseq;
if ($Tseq =~ /[Nn]{$min_gap_size}/){
    &assign_Tsubcont ($target_name, $Tseq);
}
else{
    my $count_subcon_04d = sprintf ("%04d", 1);
    my $subcon_name = $target_name . '/' . $count_subcon_04d;
    $target_subcon_seq{$subcon_name} = $Tseq;
    $target_gap_size{$subcon_name} = 0;
    $target_noscaf{$target_name} = 1;
    $total_target_subcon_len += length $Tseq;
    $target_Is{$target_name} = 1;
    $target_Ie{$target_name} = length $Tseq;
}
$Tseq = '';
close (FILE);

if (scalar keys %target_seq == scalar keys %target_noscaf){
    my $pre_ite_num = $ite_num - 1;
    print STDERR "Gaps in scaffolds have been completely closed at iteration-$pre_ite_num\n" if ($pre_ite_num > 0);
    print STDERR "Gaps in scaffolds have been already closed\n" if ($pre_ite_num == 0);
    system ("rm -rf iteration-$ite_num") if ($ite_num > 1);
    exit(1);
}


# construct subcontig information (start and end positions, and orders for each subcontig) from target subcontig position info
 
foreach my $key (keys %target_Is){
    if ($target_Is{$key} =~ /=/){
        my @start_pos = split (/=/, $target_Is{$key});
        my @end_pos = split (/=/, $target_Ie{$key});
        for (my $i = 0; $i < @start_pos; $i++){
            my $subcon_range = $start_pos[$i] . '=' . $end_pos[$i];
            push @{$Tsubcon_info{$key}}, $subcon_range;
        }
    }
    else{
        my $subcon_range = $target_Is{$key} . '=' . $target_Ie{$key};
        push @{$Tsubcon_info{$key}}, $subcon_range;
    }
}
undef %target_Is;
undef %target_Ie;

my $target_subcon_file = "$ite_dir$out_prefix-subcon.fa";

open (OUT, "> $target_subcon_file");
foreach my $name (sort keys %target_subcon_seq){
    print OUT ">$name\n";
    print OUT $target_subcon_seq{$name}, "\n";
}
close (OUT);


# get query contig sequences

my $total_query_len = 0;
open (FILE, $query_file) or die "$query_file is not found: $!\n";
while (my $line = <FILE>){
    $line =~ s/(\r\n|\n|\r)//g;
    if ($line =~ /^>/){
        if ($query_name ne '-'){
            $Qseq =~ s/^[Nn]*//;
            $Qseq =~ s/[Nn]*$//;
            $total_query_len += length $Qseq;
            $Qseq = '';
        }
        $query_name = $1 if ($line =~ /^>(\S+)/);
        die "Query name must not contain '||' or '=='" if (($query_name =~ /\|\|/) or ($query_name =~ /==/));
    }
    else{
        $Qseq .= $line;
    }
}
$Qseq =~ s/^[Nn]*//;
$Qseq =~ s/[Nn]*$//;
$total_query_len += length $Qseq;
$Qseq = '';
close (FILE);

$LR_coverage = int ($total_query_len / $target_size + 0.5) if ($LR_coverage eq 'auto');

    
# split read files

my $target_basename = basename($target_subcon_file);
my $query_basename = basename($query_file);
my $prev_prefix_base = '';
my $bowtie_path = '';
my $prev_bowtie_dir = '';
my $bowtie_dir = "$out_prefix-bowtie_align";
my $index_dir = "index";

my $align_file_flag = 0;
$align_file_flag = 1 if ((@sam_target_align_file > 0) and (@sam_query_align_file > 0));

$PE_align_tag = 1 if ($align_file_flag == 1);

if (@read_files > 0){
    if ($iteration == 1){
        if ($align_file_flag == 0){
            system ("mkdir $bowtie_dir") unless (-d $bowtie_dir);
            system ("mkdir $bowtie_dir/$index_dir") unless (-d "$bowtie_dir/$index_dir");
        }
    }
    elsif ($iteration > 1){
        system ("mkdir $ite_dir$bowtie_dir") unless (-d "$ite_dir$bowtie_dir");
        system ("mkdir $ite_dir$bowtie_dir/$index_dir") unless (-d "$ite_dir$bowtie_dir/$index_dir");
    }
}

my $total_read_num = 0;
my @read_files_2;
my @align_jobs;
my @blast_align_jobs;
my @index_jobs;

my $max_mismatch_num_coval = 1;
$max_mismatch_num_coval = int ($read_length * 0.015) if ($read_length > 100);
$max_mismatch_num_coval += 1 if ($heterozygosity == 1);


if (@read_files > 0){
    if (($read_format eq 'fastq') and ($base_call_type eq 'auto')){
        my $qual_seq;
        my $sanger_type;
        my $solexa_type;
        my $count_line = 0;
        open (FILE, $read_files[0]) or die "$read_files[0] is not found: $!\n" if ($read_files[0] !~ /\.gz$/);
        open (FILE, "gzip -dc $read_files[0] |") or die "$read_files[0] is not found: $!\n" if ($read_files[0] =~ /\.gz$/);
        while (my $line = <FILE>){
            $count_line ++;
            if ($count_line % 4 == 0){
                $qual_seq .= $line;
                last if (length $qual_seq >= 1000);
            }
        }
        $sanger_type += $qual_seq =~ /[!Ó#%&'\(\)\*,\.\/0123456789:;<=>]/g;
        $solexa_type += $qual_seq =~ /JKLMNOPQRSTUVWXYZ\[\\\]\^_`abcdefgh/g;
        close (FILE);
        if ($sanger_type > $solexa_type){
            if ($base_call_type eq 'illumina'){
                print STDERR "\n########## Warning! ##########\nYour base quality format appears a sanger type, but the job will be proceeded with illumina type.\n To change it use -qcall option.\n" if ($ite_num == 1);
            }
            else{
                $base_call_type = 'sanger';
                print STDERR "Base call quality format: sanger type\n" if ($ite_num == 1);
                print OUTLOG "Base call quality format: sanger type\n" if ($ite_num == 1);
            }
        }
        else{
            if ($base_call_type eq 'sanger'){
                print STDERR "\n########## Warning! ##########\nYour base quality format appears a illumina type, but the job will be proceeded with sanger type.\n To change use -qcall option.\n"  if ($ite_num == 1);
            }
            else{
                $base_call_type = 'illumina';
                print STDERR "Base call quality format: illumina type\n" if ($ite_num == 1);
                print OUTLOG "Base call quality format: illumina type\n" if ($ite_num == 1);
            }
        }
        $base_call_type = 'phred64' if ($base_call_type eq 'illumina');
        $base_call_type = 'phred33' if ($base_call_type eq 'sanger');
    }
    elsif ($read_format eq 'fasta'){
        my @read_file_fq;
        foreach my $file (@read_files){
            my $file_basename = basename ($file);
            my $out_fastq = "$1.mod.fastq" if ($file_basename =~ /(\S+)\.fa$/) or ($file_basename =~ /(\S+)\.fasta$/);
            $out_fastq = "$1.mod.fastq" if ($file_basename =~ /(\S+)\.fa\.gz$/) or ($file_basename =~ /(\S+)\.fasta\.gz$/);
            my $seq = '';
            my $header = '';
            open (FILE, $file) or die "$file is not found: $!\n" if ($file !~ /\.gz$/);
            open (FILE, "gzip -dc $file |") or die "$file is not found: $!\n" if ($file =~ /\.gz$/);
            open (OUT, "> $out_fastq");
            while (my $line = <FILE>){
                chomp $line;
                if ($line =~ /^>(\S+)/){
                    if ($seq ne ''){
                        print OUT "\@$header\n";
                        print OUT $seq, "\n";
                        print OUT "+\n";
                        print OUT 'H' x length ($seq), "\n";
                        $seq = '';
                    }
                    $header = $1;
                }
                else{
                    $seq .= $line;
                }
            }
            print OUT "\@$header\n";
            print OUT $seq, "\n";
            print OUT "+\n";
            print OUT 'H' x length ($seq), "\n";
            $seq = '';
            close (OUT);
            close (FILE);
            push @read_file_fq, $out_fastq;
        }
        @read_files = (@read_file_fq);
        $base_call_type = 'phred33';
    }
    if ($ite_num == 1){
        if (($thread_num / (@read_files * 6) >= 2) and ($Config{useithreads})){
            my $file_div_num = 2;
            $file_div_num = 3 if ($thread_num / (@read_files * 6) >= 3);
            $file_div_num = 4 if ($thread_num / (@read_files * 6) >= 4);
            $file_div_num *= 2 if ($ite_num > 1);
            my $wc1 = `wc -l $read_files[0]`;
            $wc1 =~ /(\d+)/;
            my $read_num1 = $1 / 4;
            my $divided_line_num = int ($read_num1 / $file_div_num) * 4;
            my $count_line = 0;
            my $count_file_num = 0;
            my @lines;
            foreach my $file (@read_files){
                open (FILE, $file) or die "$file is not found: $!\n" if ($file !~ /\.gz$/);
                open (FILE, "gzip -dc $file |") or die "$file is not found: $!\n" if ($file =~ /\.gz$/);
                while (my $line = <FILE>){
                    $count_line ++;
                    chomp $line;
                    push @lines, $line;
                    if ($count_line == $divided_line_num + 40){
                        $count_file_num ++;
                        open (OUT, "> $temp_dir/$out_prefix-read.$count_file_num.fq");
                        my $count_read = 0;
                        foreach (@lines){
                            if (($read_length >= 90) and ($insert_size < 350)){
                                $count_read ++;
                                if ($count_read == 2){
                                    my $trimmed = substr ($_, 0, 75);
                                    print OUT $trimmed, "\n";
                                }
                                elsif ($count_read == 4){
                                    my $trimmed = substr ($_, 0, 75);
                                    print OUT $trimmed, "\n";
                                    $count_read = 0;
                                }
                                else{
                                    print OUT $_, "\n";
                                }
                            }
                            else{
                                print OUT $_, "\n";
                            }
                        }
                        close (OUT);
                        push @read_files_2, "$temp_dir/$out_prefix-read.$count_file_num.fq";
                        @lines = ();
                        $count_line = 0;
                    }
                }
                if (@lines > 0){
                    $count_file_num ++;
                    open (OUT, "> $temp_dir/$out_prefix-read.$count_file_num.fq");
                    my $count_read = 0;
                    foreach (@lines){
                        if (($read_length >= 90) and ($insert_size < 350)){
                            $count_read ++;
                            if ($count_read == 2){
                                my $trimmed = substr ($_, 0, 75);
                                print OUT $trimmed, "\n";
                            }
                            elsif ($count_read == 4){
                                my $trimmed = substr ($_, 0, 75);
                                print OUT $trimmed, "\n";
                                $count_read = 0;
                            }
                            else{
                                print OUT $_, "\n";
                            }
                        }
                        else{
                            print OUT $_, "\n";
                        }
                    }
                    close (OUT);
                    push @read_files_2, "$temp_dir/$out_prefix-read.$count_file_num.fq";
                    undef @lines;
                    $count_line = 0;
                }
                close (FILE);
            }
        }
        elsif (($thread_num / (@read_files * 6) < 2) or ($Config{useithreads})){
            if (($read_length >= 90) and ($insert_size < 350)){
                my $count_file_num = 0;
                foreach my $file (@read_files){
                    $count_file_num ++;
                    open (FILE, $file) or die "$file is not found: $!\n" if ($file !~ /\.gz$/);
                    open (FILE, "gzip -dc $file |") or die "$file is not found: $!\n" if ($file =~ /\.gz$/);
                    open (OUT, "> $temp_dir/$out_prefix-read.$count_file_num.fq");
                    my $count_read = 0;
                    while (my $line = <FILE>){
                        chomp $line;
                        $count_read++;
                        if ($count_read == 2){
                            my $trimmed = substr ($line, 0, 75);
                            print OUT $trimmed, "\n";
                        }
                        elsif ($count_read == 4){
                            my $trimmed = substr ($line, 0, 75);
                            print OUT $trimmed, "\n";
                            $count_read = 0;
                        }
                        else{
                            print OUT $line, "\n";
                        }
                    }
                    close (OUT);
                    close (FILE);
                    push @read_files_2, "$temp_dir/$out_prefix-read.$count_file_num.fq";
                }
            }
            else{
                foreach my $file (@read_files){
                    if ($file =~ /(.)\.gz$/){
                        my $file_2 = $1;
                        $file_2 = basename ($file_2);
                        system ("gzip -dc $file > $file_2");
                        push @read_files_2, $file_2;
                    }
                    else{
                        push @read_files_2, $file;
                    }
                }
            }
        }
    }
    else{
        if (($thread_num / (@read_files * 6) >= 2) and ($Config{useithreads})){
            my @subread_list = <$temp_dir/$out_prefix-read.[0-9].fq>;
            @read_files_2 = (@subread_list) if (@subread_list > 0);
        }
        else{
            if (($read_length >= 90) and ($insert_size < 350)){
                my @subread_list = <$temp_dir/$out_prefix-read.[0-9].fq>;
                @read_files_2 = (@subread_list) if (@subread_list > 0);
            }
            else{
                foreach my $file (@read_files){
                    if ($file =~ /(.)\.gz$/){
                        my $file_2 = $1;
                        $file_2 = basename ($file_2);
                        system ("gzip -dc $file > $file_2");
                        push @read_files_2, $file_2;
                    }
                    else{
                        push @read_files_2, $file;
                    }
                }
            }
        }
    }
# Align between target and query contigs and align with short reads to target and query reference

    if (($thread_num > 1) and ($Config{useithreads}) and (@read_files_2 * 2 <= $thread_num)){
        my $num = 0;
        if ($ite_num == 1){
            if ((@sam_target_align_file == 0) and (@sam_query_align_file == 0)){
                my $query_index_base = "$ite_dir$bowtie_dir/$index_dir/$query_basename";
                my $thread_index1 = threads->new(\&bowtie_index, $target_subcon_file, $target_basename);
                my $thread_index2 = threads->new(\&bowtie_index, $query_file, $query_basename) if ($query_index eq '');
                push @index_jobs, $thread_index1, $thread_index2 if ($query_index eq '');
                push @index_jobs, $thread_index1 if ($query_index ne '');
                foreach (@index_jobs){
                    $_->join;
                }
                undef @index_jobs;
                foreach my $reads (@read_files_2){
                    $num ++;
                    my ($thread_t) = threads->new(\&read_align, $reads, $num, 'T');
                    my ($thread_q) = threads->new(\&read_align, $reads, $num, 'Q');
                    push @align_jobs, $thread_t, $thread_q;
                }
            }
            elsif (@sam_query_align_file == 0){
                my $thread_index1 = threads->new(\&bowtie_index, $query_file, $query_basename);
                push @index_jobs, $thread_index1;
                foreach (@index_jobs){
                    $_->join;
                }
                undef @index_jobs;
                foreach my $reads (@read_files_2){
                    $num ++;
                    my ($thread_q) = threads->new(\&read_align, $reads, $num, 'Q');
                    push @align_jobs, $thread_q;
                }
            }
            elsif (@sam_target_align_file == 0){
                my $thread_index1 = threads->new(\&bowtie_index, $target_subcon_file, $target_basename);
                push @index_jobs, $thread_index1;
                foreach (@index_jobs){
                    $_->join;
                }
                undef @index_jobs;
                foreach my $reads (@read_files_2){
                    $num ++;
                    my ($thread_t) = threads->new(\&read_align, $reads, $num, 'T');
                    push @align_jobs, $thread_t;
                }
            }
        }
        else{
            my $thread_index1 = threads->new(\&bowtie_index, $target_subcon_file, $target_basename);
            push @index_jobs, $thread_index1;
            foreach (@index_jobs){
                $_->join;
            }
            undef @index_jobs;
            foreach my $reads (@read_files_2){
                $num ++;
                my ($thread_t) = threads->new(\&read_align, $reads, $num, 'T');
                push @align_jobs, $thread_t;
            }
        }
        undef @read_files_2;
        foreach (@align_jobs){
            my ($num1, $file_name, $tag) = $_->join;
            $align_file = $file_name if ($tag eq 'N');
            push @sam_target_align_file, $file_name if ($tag eq 'T');
            push @sam_query_align_file, $file_name if ($tag eq 'Q');
            print STDERR "Bowtie alignment of read file $num1 with target contigs was completed:\n" if ($num1 > 0) and ($tag eq 'T');
            print STDERR "Bowtie alignment of read file $num1 with query contigs was completed:\n" if ($num1 > 0) and ($tag eq 'Q');
        }
        undef @align_jobs;
    }
    else{
        print STDERR "Your Perl has not been compiled to use 'threads' module. PE-read alignments will be processed at a single CPU.\n" if (!$Config{useithreads} and $thread_num > 1);
        my $num = 0;
        if ($ite_num == 1){
            if ((@sam_target_align_file == 0) and (@sam_query_align_file == 0)){
                system ("bowtie2-build -q -f $target_subcon_file $ite_dir$bowtie_dir/$index_dir/$target_basename 2>/dev/null");
                system ("bowtie2-build -q -f $query_file $ite_dir$bowtie_dir/$index_dir/$query_basename 2>/dev/null");
                foreach my $read_file (@read_files){
                    $num ++;
                    my $read_file_2 = $read_file;
                    $read_file_2 = basename ($read_file_2);
                    if ($read_file =~ /(.+)\.gz$/){
                        $read_file_2 = $1;
                        system ("gzip -dc $read_file > $read_file_2");
                    }
                    system ("bowtie2 --$base_call_type -p $thread_num -x $ite_dir$bowtie_dir/$index_dir/$target_basename $read_file_2 | $Bin/coval-filter-short.pl -n $max_mismatch_num_coval -r $target_subcon_file - > $ite_dir$bowtie_dir/$out_prefix-target-align-se-$num.sam");
                    print STDERR "Bowtie alignment of all read files with target was finished:\n";
                    push @sam_target_align_file, "$ite_dir$bowtie_dir/$out_prefix-target-align-se-$num.sam";
                    my $index_base = "$ite_dir$bowtie_dir/$index_dir/$query_basename";
                    $index_base = $query_index if ($query_index ne '');
                    system ("bowtie2 --$base_call_type -k $LR_coverage -p $thread_num -x $index_base $read_file_2 | $Bin/coval-filter-short.pl -n $max_mismatch_num_coval -r $query_file - > $ite_dir$bowtie_dir/$out_prefix-query-align-se-$num.sam");
                    print STDERR "Bowtie alignment of all read files with query was finished:\n";
                    push @sam_query_align_file, "$ite_dir$bowtie_dir/$out_prefix-query-align-se-$num.sam";
                }
            }
            elsif (@sam_query_align_file == 0){
                system ("bowtie2-build -q -f $query_file $ite_dir$bowtie_dir/$index_dir/$query_basename 2>/dev/null");
                foreach my $read_file (@read_files){
                    $num ++;
                    my $read_file_2 = $read_file;
                    $read_file_2 = basename ($read_file_2);
                    if ($read_file =~ /(.+)\.gz$/){
                        $read_file_2 = $1;
                        system ("gzip -dc $read_file > $read_file_2");
                    }
                    my $index_base = "$ite_dir$bowtie_dir/$index_dir/$query_basename";
                    $index_base = $query_index if ($query_index ne '');
                    system ("bowtie2 --$base_call_type -k $LR_coverage -p $thread_num -x $index_base $read_file_2 | coval-filter-short.pl -n $max_mismatch_num_coval -r $query_file - > $ite_dir$bowtie_dir/$out_prefix-query-align-se-$num.sam");
                    print STDERR "Bowtie alignment of all read files with query was finished:\n";
                    push @sam_query_align_file, "$ite_dir$bowtie_dir/$out_prefix-query-align-se-$num.sam";
                }
            }
            elsif (@sam_target_align_file == 0){
                system ("bowtie2-build -q -f $target_subcon_file $ite_dir$bowtie_dir/$index_dir/$target_basename 2>/dev/null");
                foreach my $read_file (@read_files){
                    $num ++;
                    my $read_file_2 = $read_file;
                    $read_file_2 = basename ($read_file_2);
                    if ($read_file =~ /(.+)\.gz$/){
                        $read_file_2 = $1;
                        system ("gzip -dc $read_file > $read_file_2");
                    }
                    system ("bowtie2 --$base_call_type -p $thread_num -x $ite_dir$bowtie_dir/$index_dir/$target_basename $read_file_2 | $Bin/coval-filter-short.pl -n $max_mismatch_num_coval -r $target_subcon_file - > $ite_dir$bowtie_dir/$out_prefix-target-align-se-$num.sam");
                    print STDERR "Bowtie alignment of all read files with target was finished:\n";
                    push @sam_target_align_file, "$ite_dir$bowtie_dir/$out_prefix-target-align-se-$num.sam";
                }
            }
        }
        else{
            system ("bowtie2-build -q -f $target_subcon_file $ite_dir$bowtie_dir/$index_dir/$target_basename 2>/dev/null");
            foreach my $read_file (@read_files){
                $num ++;
                my $read_file_2 = $read_file;
                if ($read_file =~ /(.+)\.gz$/){
                    $read_file_2 = $1;
                    $read_file_2 = basename ($read_file_2);
                    system ("gzip -dc $read_file > $read_file_2");
                }
                system ("bowtie2 --$base_call_type -p $thread_num -x $ite_dir$bowtie_dir/$index_dir/$target_basename $read_file_2 | $Bin/coval-filter-short.pl -n $max_mismatch_num_coval -r $target_subcon_file - > $ite_dir$bowtie_dir/$out_prefix-target-align-se-$num.sam");
                print STDERR "Bowtie alignment of all read files with target was finished:\n";
                push @sam_target_align_file, "$ite_dir$bowtie_dir/$out_prefix-target-align-se-$num.sam";
            }
        }
    }
#    my @file_list = <$temp_dir/$out_prefix-read.[0-9].fq>;
#        system ("rm -f $temp_dir/$out_prefix-read.[0-9].fq") if (@file_list > 0);
}

if ($ite_num > 1){
    my @query_sam_list = <iteration-1/$bowtie_dir/$out_prefix-query-align-se-[0-9].sam>;
    @sam_query_align_file = (@query_sam_list) if (@query_sam_list > 0);
}

my %query_seq;
$query_name = '';
$Qseq = '';
open (FILE, $query_file) or die "$query_file is not found: $!\n";
while (my $line = <FILE>){
    $line =~ s/(\r\n|\n|\r)//g;
    if ($line =~ /^>/){
        if ($query_name ne '-'){
            $Qseq =~ s/^[Nn]*//;
            $Qseq =~ s/[Nn]*$//;
            $query_seq{$query_name} = uc $Qseq;
            $Qseq = '';
        }
        $query_name = $1 if ($line =~ /^>(\S+)/);
        die "Query name must not contain '||' or '=='" if (($query_name =~ /\|\|/) or ($query_name =~ /==/));
    }
    else{
        $Qseq .= $line;
    }
}
$Qseq =~ s/^[Nn]*//;
$Qseq =~ s/[Nn]*$//;
$query_seq{$query_name} = uc $Qseq;
$Qseq = '';
close (FILE);

if (($align_file eq '') or ($align_file_Q eq '')){
    if ($ite_num == 1){
        if ((($align_file eq '') and ($align_file_Q eq '')) or (($align_file eq '') and ($iteration > 1))){
            system ("mkdir $blast_dir") unless (-d $blast_dir);
            system ("cp -f $target_subcon_file $blast_dir/");
            system ("makeblastdb -in $blast_dir/$target_basename -dbtype nucl -logfile $temp_dir/blastdb.log");
            system ("cp -f $query_file $blast_dir/");
            system ("makeblastdb -in $blast_dir/$query_basename -dbtype nucl -logfile $temp_dir/blastdb.log");
        }
        elsif ($align_file eq ''){
            system ("mkdir $blast_dir") unless (-d $blast_dir);
            system ("cp -f $target_subcon_file $blast_dir/");
            system ("makeblastdb -in $blast_dir/$target_basename -dbtype nucl -logfile $temp_dir/blastdb.log");
        }
        elsif ($align_file_Q eq ''){
            system ("mkdir $blast_dir") unless (-d $blast_dir);
            system ("cp -f $query_file $blast_dir/");
            system ("makeblastdb -in $blast_dir/$query_basename -dbtype nucl -logfile $temp_dir/blastdb.log");
        }
    }
    else{
        system ("mkdir $blast_dir") unless (-d $blast_dir);
        system ("cp -f $target_subcon_file $blast_dir/");
        system ("makeblastdb -in $blast_dir/$target_basename -dbtype nucl -logfile $temp_dir/blastdb.log");
    }
    if ($thread_num >= 2 and $Config{useithreads}){
        if ((($thread_num >= 4) and ($ite_num == 1)) or (($thread_num >= 2) and ($ite_num > 1))){
            if ((($thread_num >= 8) and ($ite_num == 1)) or (($thread_num >= 4) and ($ite_num > 1))){
                my $query_num = scalar keys %query_seq;
                my $div_num = int ($thread_num / 4 + 0.5);
                $div_num = int ($thread_num / 2 + 0.5) if ($ite_num > 1) or (($ite_num == 1) and (($align_file ne '') or ($align_file_Q ne '')));
                $div_num = 24 if ($div_num > 24);
                my $div_lines = int ($query_num / $div_num) + 1;
                my $count_div = 1;
                my $count_q = 0;
                my %query_subseq;
                my @query_subfiles;
                foreach my $qname (keys %query_seq){
                    $count_q ++;
                    $query_subseq{$qname} = $query_seq{$qname};
                    if ($count_q == $div_lines){
                        open (OUT, ">$temp_dir/$out_prefix-query-$count_div.fa");
                        foreach my $qname2 (keys %query_subseq){
                            print OUT ">$qname2\n";
                            print OUT $query_subseq{$qname2}, "\n";
                        }
                        close (OUT);
                        push @query_subfiles, "$temp_dir/$out_prefix-query-$count_div.fa";
                        $count_div ++;
                        %query_subseq = ();
                        $count_q = 0;
                    }
                }
                if ($count_q > 0){
                    open (OUT, ">$temp_dir/$out_prefix-query-$count_div.fa");
                    foreach my $qname2 (keys %query_subseq){
                        print OUT ">$qname2\n";
                        print OUT $query_subseq{$qname2}, "\n";
                    }
                    close (OUT);
                    push @query_subfiles, "$temp_dir/$out_prefix-query-$count_div.fa";
                    undef %query_subseq;
                }
                undef %query_seq;
                my $subfile_num = 0;
                foreach (@query_subfiles){
                    $subfile_num ++;
                    my $align_Tsub_out = "$temp_dir/$out_prefix-T-$subfile_num.blast.out";
                    my $align_Qsub_out = "$temp_dir/$out_prefix-Q-$subfile_num.blast.out";
                    if ($ite_num == 1){
                        if (($align_file eq '') and ($align_file_Q eq '')){
                            my ($thread_t) = threads->new(\&blast_align_thread, $_, $align_Tsub_out);
                            my ($thread_q) = threads->new(\&blast_align_thread_Q, $_, $align_Qsub_out);
                            push @blast_align_jobs, $thread_t, $thread_q;
                        }
                        elsif ($align_file eq ''){
                            my ($thread_t) = threads->new(\&blast_align_thread, $_, $align_Tsub_out);
                            push @blast_align_jobs, $thread_t;
                        }
                        elsif ($align_file_Q eq ''){
                            my ($thread_q) = threads->new(\&blast_align_thread_Q, $_, $align_Qsub_out);
                            push @blast_align_jobs, $thread_q;
                        }
                    }
                    else{
                        my ($thread_t) = threads->new(\&blast_align_thread, $_, $align_Tsub_out);
                        push @blast_align_jobs, $thread_t;
                    }
                }
                undef @query_subfiles;
                foreach (@blast_align_jobs){
                    my ($file_name, $tag) = $_->join;
                    push @blast_target_align_file, $file_name if ($tag eq 'T');
                    push @blast_query_align_file, $file_name if ($tag eq 'Q');
                    print STDERR "BLASTn alignmnet for $file_name --- completed:\n";
                }
                if ($align_file eq ''){
                    my $align_arg = '';
                    $align_file = "$ite_dir$out_prefix-1m.blast.out";
                    foreach (@blast_target_align_file){
                        $align_arg .= $_ . ' ';
                    }
                    system ("cat $align_arg > $align_file");
                }
                if (($ite_num == 1) and ($align_file_Q eq '')){
                    my $align_arg_Q = '';             
                    $align_file_Q = "$ite_dir$out_prefix-QQ.blast.out";
                    foreach (@blast_query_align_file){
                        $align_arg_Q .= $_ . ' ';
                    }
                    system ("cat $align_arg_Q > $align_file_Q");
                }
            }
            else{
                if (($align_file eq '') and ($align_file_Q eq '')){
                    $align_file = "$ite_dir$out_prefix-1m.blast.out";
                    $align_file_Q = "$ite_dir$out_prefix-QQ.blast.out";
                    if ($ite_num == 1){
                        my ($thread_t) = threads->new(\&blast_align_thread, $query_file, $align_file);
                        my ($thread_q) = threads->new(\&blast_align_thread_Q, $query_file, $align_file_Q);
                        push @blast_align_jobs, $thread_t, $thread_q;
                    }
                    else{
                        my ($thread_t) = threads->new(\&blast_align_thread, $query_file, $align_file);
                        push @blast_align_jobs, $thread_t;
                    }
                }
                elsif ($align_file eq ''){
                    $align_file = "$ite_dir$out_prefix-1m.blast.out";
                    my ($thread_t) = threads->new(\&blast_align_thread, $query_file, $align_file);
                    push @blast_align_jobs, $thread_t;
                }
                elsif (($align_file_Q eq '') and ($ite_num == 1)){
                    $align_file_Q = "$ite_dir$out_prefix-QQ.blast.out";
                    my ($thread_q) = threads->new(\&blast_align_thread_Q, $query_file, $align_file_Q);
                    push @blast_align_jobs, $thread_q;
                }
                foreach (@blast_align_jobs){
                    my ($file_name, $tag) = $_->join;
                    print STDERR "BLASTn alignmnet for target-query --- completed:\n" if ($tag eq 'T');
                    print STDERR "BLASTn alignmnet for query-query --- completed:\n" if ($tag eq 'Q');
                }
            }
            undef @blast_align_jobs;
            my @query_subfile_list = <$temp_dir/$out_prefix-query-[0-9].fa>;
            my @blast_subfile_list = <$temp_dir/$out_prefix-[TQ]-[0-9].blast.out>;
            system ("rm -f $temp_dir/$out_prefix-query-[0-9].fa") if (@query_subfile_list > 0);
            system ("rm -f $temp_dir/$out_prefix-[TQ]-[0-9].blast.out") if (@blast_subfile_list > 0);
        }
        else{
            $align_file = &blast_align () if ($align_file eq '');
            $align_file_Q = &blast_align_Q () if ($ite_num == 1) and ($align_file_Q eq '');
        }
    }
    else{
        $align_file = &blast_align () if ($align_file eq '');
        $align_file_Q = &blast_align_Q () if ($ite_num == 1) and ($align_file_Q eq '');
    }
}


sub blast_align_thread{
    my $query = shift;
    my $align_out = shift;
    system ("blastn -db $blast_dir/$target_basename -query $query -out $align_out -outfmt 6 -num_alignments 100 -perc_identity 95 -num_threads 2");
    threads->yield();
    sleep 1;
    return ($align_out, 'T');
}

sub blast_align_thread_Q{
    my $query = shift;
    my $align_out = shift;
    system ("blastn -db $blast_dir/$query_basename -query $query -out $align_out -outfmt 6 -num_alignments 100 -perc_identity 95 -num_threads 2");
    threads->yield();
    sleep 1;
    return ($align_out, 'Q');
}

sub blast_align{
    my $align_file = "$ite_dir$out_prefix-1m.blast.out";
    my $blast_dir = 'db';
    system ("mkdir $blast_dir") unless (-d $blast_dir);
    system ("cp -f $target_subcon_file $blast_dir/");
    system ("makeblastdb -in $blast_dir/$target_basename -dbtype nucl -logfile $temp_dir/blastdb.log");
    system ("blastn -db $blast_dir/$target_basename -query $query_file -out $align_file -outfmt 6 -num_alignments 100 -perc_identity 95 -num_threads 2") if ($Config{useithreads});
    system ("blastn -db $blast_dir/$target_basename -query $query_file -out $align_file -outfmt 6 -num_alignments 100 -perc_identity 95") unless ($Config{useithreads});
    print STDERR "BLAST alignment for target-query --- completed:\n";
    return ($align_file);
}

sub blast_align_Q{
    my $align_file = "$ite_dir$out_prefix-QQ.blast.out";
    my $blast_dir = 'db';
    system ("mkdir $blast_dir") unless (-d $blast_dir);
    system ("cp -f $query_file $blast_dir/");
    system ("makeblastdb -in $blast_dir/$query_basename -dbtype nucl -logfile $temp_dir/blastdb.log");
    system ("blastn -db $blast_dir/$query_basename -query $query_file -out $align_file -outfmt 6 -num_alignments 100 -perc_identity 95 -num_threads 2") if ($Config{useithreads});
    system ("blastn -db $blast_dir/$query_basename -query $query_file -out $align_file -outfmt 6 -num_alignments 100 -perc_identity 95") unless ($Config{useithreads});
    print STDERR "BLAST alignment for query-query --- completed:\n";
    return ($align_file);
}

sub bowtie_index{
    my $ref_file = shift;
    my $ref_base = shift;
    system ("bowtie2-build -q -f $ref_file $ite_dir$bowtie_dir/$index_dir/$ref_base 2>/dev/null");
    threads->yield();
    sleep 1;
}

sub read_align{
    my $readfile = shift;
    my $readnum = shift;
    my $ref_tag = shift;
    if ($ref_tag eq 'T'){
        system ("bowtie2 --$base_call_type -x $ite_dir$bowtie_dir/$index_dir/$target_basename -p 3 $readfile | $Bin/coval-filter-short.pl -n $max_mismatch_num_coval -r $target_subcon_file - > $ite_dir$bowtie_dir/$out_prefix-target-align-se-$readnum.sam") if ($ite_num == 1);
        system ("bowtie2 --$base_call_type -x $ite_dir$bowtie_dir/$index_dir/$target_basename -p 5 $readfile | $Bin/coval-filter-short.pl -n $max_mismatch_num_coval -r $target_subcon_file - > $ite_dir$bowtie_dir/$out_prefix-target-align-se-$readnum.sam") if ($ite_num > 1);
    }
    elsif ($ref_tag eq 'Q'){
        my $index_base = "$ite_dir$bowtie_dir/$index_dir/$query_basename";
        $index_base = $query_index if ($query_index ne '');
        system ("bowtie2 --$base_call_type -k $LR_coverage -x $index_base -p 3 $readfile | $Bin/coval-filter-short.pl -n $max_mismatch_num_coval -r $query_file - > $ite_dir$bowtie_dir/$out_prefix-query-align-se-$readnum.sam");
    }
    threads->yield();
    sleep 1;
    return ($readnum, "$ite_dir$bowtie_dir/$out_prefix-target-align-se-$readnum.sam", 'T') if ($ref_tag eq 'T');
    return ($readnum, "$ite_dir$bowtie_dir/$out_prefix-query-align-se-$readnum.sam", 'Q') if ($ref_tag eq 'Q');
}

if (scalar keys %query_seq == 0){
    $query_name = '';
    $Qseq = '';
    open (FILE, $query_file) or die "$query_file is not found: $!\n";
    while (my $line = <FILE>){
        $line =~ s/(\r\n|\n|\r)//g;
        if ($line =~ /^>/){
            if ($query_name ne '-'){
                $Qseq =~ s/^[Nn]*//;
                $Qseq =~ s/[Nn]*$//;
                $query_seq{$query_name} = uc $Qseq;
                $Qseq = '';
            }
            $query_name = $1 if ($line =~ /^>(\S+)/);
            die "Query name must not contain '||' or '=='" if (($query_name =~ /\|\|/) or ($query_name =~ /==/));
        }
        else{
            $Qseq .= $line;
        }
    }
    $Qseq =~ s/^[Nn]*//;
    $Qseq =~ s/[Nn]*$//;
    $query_seq{$query_name} = uc $Qseq;
    $Qseq = '';
    close (FILE);
}


my %T_match_5;     # query sequences that cover the 5'-end regions of target subsequences (subcontigs)
my %T_match_3;     # query sequences that cover the 3'-end regions of target subsequences
my %T_match_A;     # query sequences that cover all the region of target subsequences
my %T_match_I;     # query sequences that are covered within target subsequences
my %T_match_ALL;   # query sequences that cover target subsequences (%T_match_5 + %T_match_3 + %T_match_A)
my %TQ_cover;      # target subseqeunces in the original target scaffolds that are covered with queries
my %Q_match;       # total length of unmatch regions of a query aligned with target seqeunces
my %TQ_strand;     # strand of queries that cover target subseqeunces
my %included_qname;  # query name merged within target seqeunces
my %included_tname;  # target name merged within query seqeunces
my %included_qname_Q;
my %included_tname_Q;
my %Q_match_5;
my %Q_match_3;
my %evaluate_pairs;
my $sum_score = 0;
my $num_score = 0;
my $sum_score_2 = 0;
my $num_score_2 = 0;
my %target_info;

my %target_read_pos_F;                               # store alignment positions for each read from sam files for target scaffolds
my %target_read_pos_R;
my %query_read_pos_F;
my %query_read_pos_R;
my $count = 0;

my $max_insert_size = $insert_size + $SD_insert * 5;
my $min_insert_size = $insert_size - $SD_insert * 5;

my $target_read_num = 0;
my $query_read_num = 0;

if (@sam_target_align_file > 0){            # store alignment positions for each read from sam files for target contigs
    foreach my $sam_file (@sam_target_align_file){
        &get_read_pos ($sam_file, 'T');
    }
}
if ((@sam_query_align_file > 0) and ($ite_num == 1)){
    foreach my $sam_file (@sam_query_align_file){
        &get_read_pos ($sam_file, 'Q');
    }
}

if ($ite_num > 1){
    open (FILE, "Query-Query-align/$out_prefix-query-F.txt") or die "Query-Query-align/$out_prefix-query-F.txt is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my @line = split (/\t/, $line);
        ${$query_read_pos_F{$line[0]}}{$line[1]} = $line[2];
    }
    close (FILE);
    open (FILE, "Query-Query-align/$out_prefix-query-R.txt") or die "Query-Query-align/$out_prefix-query-R.txt is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my @line = split (/\t/, $line);
        ${$query_read_pos_R{$line[0]}}{$line[1]} = $line[2];
    }
    close (FILE);
}

my $score_file = "$ite_dir$out_prefix-score.txt";
if ($PE_align_tag == 1){
    my $ave_target_coverage = int ($target_read_num * $read_length / $total_target_subcon_len * 10) / 10;
    print STDERR "Total number of target subcontig sequences: ", scalar keys %target_subcon_seq, "\n";
    print STDERR "Apparent alignment coverage for target: $ave_target_coverage x\n";
    print OUTLOG "Total number of target subcontig sequences: ", scalar keys %target_subcon_seq, "\n";
    print OUTLOG "Apparent alignment coverage for target: $ave_target_coverage x\n";
    my $ave_query_coverage = 0;
    if ($ite_num == 1){
        $ave_query_coverage = int ($query_read_num * $read_length / $total_query_len * 10) / 10;
        print STDERR "Total number of query sequences: ", scalar keys %query_seq, "\n";
        print STDERR "Apparent alignment coverage for query: $ave_query_coverage x\n";
        print OUTLOG "Total number of query sequences: ", scalar keys %query_seq, "\n";
        print OUTLOG "Apparent alignment coverage for query: $ave_query_coverage x\n";
    }
    $min_insert_size = $read_length * 2 if ($min_insert_size < $read_length * 2);
    open OUT, "> $score_file";
    print OUT "##average-target-coverage: $ave_target_coverage\t##average-query-coverage: $ave_query_coverage\n" if ($ite_num == 1);
    print OUT "##average-target-coverage: $ave_target_coverage\n" if ($ite_num > 1);
    close (OUT);
    close (FILE);
}

if (($iteration > 1) and ($ite_num == 1)){
    system ("mkdir Query-Query-align") unless (-d "Query-Query-align");
    my $out_query_aligned_read_F = "Query-Query-align/$out_prefix-query-F.txt";
    my $out_query_aligned_read_R = "Query-Query-align/$out_prefix-query-R.txt";
    open (OUT, "> $out_query_aligned_read_F");
    foreach my $key1 (keys %query_read_pos_F){
        foreach my $key2 (keys %{$query_read_pos_F{$key1}}){
            print OUT "$key1\t$key2\t${$query_read_pos_F{$key1}}{$key2}\n";
        }
    }
    close (OUT);
    open (OUT, "> $out_query_aligned_read_R");
    foreach my $key1 (keys %query_read_pos_R){
        foreach my $key2 (keys %{$query_read_pos_R{$key1}}){
            print OUT "$key1\t$key2\t${$query_read_pos_R{$key1}}{$key2}\n";
        }
    }
    close (OUT);
}


sub get_read_pos{
    my $file = shift;
    my $tag = shift;
    open (FILE, $file) or die "$file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^\@/) or ($line =~ /^SSAHA/) or ($line =~ /^$/);
        my @line = split(/\s+/, $line);
        next if (@line < 6);
        next if ($line[2] eq '*');
        my $ref_name = $line[2];
        my $ref_pos = $line[3];
        my $read_name = $line[0];
        my $read_dir = '';
        if (($line[1] == 0) or ($line[1] == 256)){
            $read_dir = 'F';
        }
        elsif (($line[1] == 16) or ($line[1] == 272)){
            $read_dir = 'R';
        }
        else{
            $read_dir = 'FR';
        }
        if ($tag eq 'T'){
            $target_read_num ++;
            if ($read_dir eq 'F'){
                ${$target_read_pos_F{$ref_name}}{$ref_pos} .= $read_name . '=';
            }
            elsif ($read_dir eq 'R'){
                ${$target_read_pos_R{$ref_name}}{$ref_pos} .= $read_name . '=';
            }
            else{
                ${$target_read_pos_F{$ref_name}}{$ref_pos} .= $read_name . '=';
                ${$target_read_pos_R{$ref_name}}{$ref_pos} .= $read_name . '=';
            }
        }
        elsif ($tag eq 'Q'){
            $query_read_num ++;
            if ($read_dir eq 'F'){
                ${$query_read_pos_F{$ref_name}}{$ref_pos} .= $read_name . '=';
            }
            elsif ($read_dir eq 'R'){
                ${$query_read_pos_R{$ref_name}}{$ref_pos} .= $read_name . '=';
            }
            else{
                ${$query_read_pos_F{$ref_name}}{$ref_pos} .= $read_name . '=';
                ${$query_read_pos_R{$ref_name}}{$ref_pos} .= $read_name . '=';
            }
        }
    }
    close (FILE);
}


# get consistent alignment info from a blast alignment file

$target_name = '';
$query_name = '';

my $out_blast_filt = "$ite_dir$out_prefix-2m.blast-filt.txt";
open (OUT1, "> $out_blast_filt");
close (OUT1);

&blast_align_assign ($align_file, 'T');

my $blast_QQ5_out = "Query-Query-align/$out_prefix-QQ5.align.txt";
my $blast_QQ3_out = "Query-Query-align/$out_prefix-QQ3.align.txt";
my $blast_QQincQ_out = "Query-Query-align/$out_prefix-QQincQ.align.txt";
my $blast_QQincT_out = "Query-Query-align/$out_prefix-QQincT.align.txt";

if ($iteration > 1){
    if ($ite_num == 1){
        system ("mkdir Query-Query-align") unless (-d "Query-Query-align");
        open (OUT1, "> $ite_dir$out_prefix-query-query.blast-filt.txt");
        close (OUT1);
        &blast_align_assign ($align_file_Q, 'Q');
        
        open (OUT, "> $blast_QQ5_out");
        foreach my $key1 (keys %Q_match_5){
            foreach my $key2 (keys %{$Q_match_5{$key1}}){
                print OUT "$key1\t$key2\t${$Q_match_5{$key1}}{$key2}\n";
            }
        }
        close (OUT);
        open (OUT, "> $blast_QQ3_out");
        foreach my $key1 (keys %Q_match_3){
            foreach my $key2 (keys %{$Q_match_3{$key1}}){
                print OUT "$key1\t$key2\t${$Q_match_3{$key1}}{$key2}\n";
            }
        }
        close (OUT);
        open (OUT, "> $blast_QQincQ_out");
        foreach my $key1 (keys %included_qname_Q){
            foreach my $key2 (keys %{$included_qname_Q{$key1}}){
                print OUT "$key1\t$key2\t${$included_qname_Q{$key1}}{$key2}\n";
            }
        }
        close (OUT);
        open (OUT, "> $blast_QQincT_out");
        foreach my $key1 (keys %included_tname_Q){
            foreach my $key2 (keys %{$included_tname_Q{$key1}}){
                print OUT "$key1\t$key2\t${$included_tname_Q{$key1}}{$key2}\n";
            }
        }
        close (OUT);
    }
    else{
        open (FILE, $blast_QQ5_out) or die "$blast_QQ5_out is not found: $!\n";
        while (my $line = <FILE>){
            chomp $line;
            my @line = split (/\t/, $line);
            ${$Q_match_5{$line[0]}}{$line[1]} = $line[2];
        }
        close (FILE);
        open (FILE, $blast_QQ3_out) or die "$blast_QQ3_out is not found: $!\n";
        while (my $line = <FILE>){
            chomp $line;
            my @line = split (/\t/, $line);
            ${$Q_match_3{$line[0]}}{$line[1]} = $line[2];
        }
        close (FILE);
        open (FILE, $blast_QQincQ_out) or die "$blast_QQincQ_out is not found: $!\n";
        while (my $line = <FILE>){
            chomp $line;
            my @line = split (/\t/, $line);
            ${$included_qname_Q{$line[0]}}{$line[1]} = $line[2];
        }
        close (FILE);
        open (FILE, $blast_QQincT_out) or die "$blast_QQincT_out is not found: $!\n";
        while (my $line = <FILE>){
            chomp $line;
            my @line = split (/\t/, $line);
            ${$included_tname_Q{$line[0]}}{$line[1]} = $line[2];
        }
        close (FILE);
    }
}

sub blast_align_assign {
    my ($align_file_2, $tag) = @_;
    my %partial_match_5;
    open (FILE, $align_file_2) or die "$align_file_2 is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my @line = split (/\s+/, $line);
        my $qname = $line[0];
        my $tname = $line[1];
        my $identity = $line[2];
        my $overlen = $line[3];
        my $qstart = $line[6];
        my $qend = $line[7];
        my $strand = 'Plus' if ($line[8] <= $line[9]);
        $strand = 'Minus' if ($line[8] > $line[9]);
        my $tstart = $line[8] if ($strand eq 'Plus');
        $tstart = $line[9] if ($strand eq 'Minus');
        my $tend = $line[9] if ($strand eq 'Plus');
        $tend = $line[8] if ($strand eq 'Minus');
        my $qlen = length $query_seq{$qname};
        my $tlen = length $target_subcon_seq{$tname} if ($tag eq 'T');
        $tlen = length $query_seq{$tname} if ($tag eq 'Q');
        my $qrange = "$qstart=$qend";
        if ((($identity >= 97) and ($overlen >= 50) and ($heterozygosity == 0) and ($qname ne $tname)) or (($identity >= 96) and ($overlen >= 50) and ($heterozygosity == 1) and ($qname ne $tname))){
            if ((($qstart <= 3) and ($qend >= $qlen - 2)) or (($tstart <= 3) and ($tend >= $tlen - 2))){
                my $align_info = "$qname==$qstart==$qend==$tname==$tstart==$tend==$strand==$overlen==$identity";
                if (!exists ${$target_info{$tname}}{$qname}){
                    ${$target_info{$tname}}{$qname} = $align_info;
                }
            }
            elsif ((($qstart <= 3) and ($tend >=  $tlen - 2) and ($strand eq 'Plus')) or (($tstart <= 3) and ($qend >= $qlen - 2) and ($strand eq 'Plus')) or (($qstart <= 3) and ($tstart <= 3) and ($strand eq 'Minus')) or (($qend >= $qlen - 2) and ($tend >=  $tlen - 2) and ($strand eq 'Minus'))){
                my $align_info = "$qname==$qstart==$qend==$tname==$tstart==$tend==$strand==$overlen==$identity";
                if (!exists ${$target_info{$tname}}{$qname}){
                    ${$target_info{$tname}}{$qname} = $align_info;
                }
            }
            if (!exists ${$target_info{$tname}}{$qname}){
                ${${${$partial_match_5{$tname}}{$qname}}{$strand}}{$tstart} = "$tend==$qstart==$qend==$overlen==$identity";
            }
        }
    }
    close (FILE);
    
    foreach my $tname (keys %partial_match_5){
        my $tlen = length $target_subcon_seq{$tname} if ($tag eq 'T');
        $tlen = length $query_seq{$tname} if ($tag eq 'Q');
        foreach my $qname (keys %{$partial_match_5{$tname}}){
            my $qlen = length $query_seq{$qname};
            foreach my $strand (keys %{${$partial_match_5{$tname}}{$qname}}){
                my $pre_tend = 0;
                my $pre_tstart = 0;
                my $pre_qend = 0;
                my $pre_qstart = 0;
                my $pre_overlen = 0;
                my $pre_ident = 0;
                foreach my $tstart (sort {$a<=>$b} keys %{${${$partial_match_5{$tname}}{$qname}}{$strand}}){
                    my ($tend, $qstart, $qend, $overlen, $ident) = split (/==/, ${${${$partial_match_5{$tname}}{$qname}}{$strand}}{$tstart});
                    if ($pre_tend > 0){
                        my $t_distance = abs ($tstart - $pre_tend);
                        my $q_distance = abs ($qstart - $pre_qend) if ($strand eq 'Plus');
                        $q_distance = abs ($qend - $pre_qstart) if ($strand eq 'Minus');
                        if (($q_distance <= $max_indel_size) and ($t_distance <= $max_indel_size)){
                            my $new_ident = int (($pre_ident * $pre_overlen + $ident * $overlen) / ($pre_overlen + $overlen));
                            my $new_overlen = $tend - $pre_tstart + 1;
                            $pre_tend = $tend;
                            $pre_qend = $qend if ($strand eq 'Plus');
                            $pre_qstart = $qstart if ($strand eq 'Minus');
                            $pre_ident = $new_ident;
                            $pre_overlen = $new_overlen;
                            next;
                        }
                        else{
                            if ((($pre_qstart > 0) and ($pre_qstart <= 3) and ($pre_qend >= $qlen - 2)) or (($pre_tstart > 0) and ($pre_tstart <= 3) and ($pre_tend >= $tlen - 2))){
                                my $align_info = "$qname==$pre_qstart==$pre_qend==$tname==$pre_tstart==$pre_tend==$strand==$pre_overlen==$pre_ident";
                                if (!exists ${$target_info{$tname}}{$qname}){
                                    ${$target_info{$tname}}{$qname} = $align_info;
                                }
                            }
                            elsif ((($pre_qstart > 0) and ($pre_qstart <= 3) and ($pre_tend >=  $tlen - 2) and ($strand eq 'Plus')) or (($pre_tstart > 0) and ($pre_tstart <= 3) and ($pre_qend >= $qlen - 2) and ($strand eq 'Plus')) or (($pre_qstart > 0) and ($pre_qstart <= 3) and ($pre_tstart <= 3) and ($strand eq 'Minus')) or (($pre_qend >= $qlen - 2) and ($pre_tend >=  $tlen - 2) and ($strand eq 'Minus'))){
                                my $align_info = "$qname==$pre_qstart==$pre_qend==$tname==$pre_tstart==$pre_tend==$strand==$pre_overlen==$pre_ident";
                                if (!exists ${$target_info{$tname}}{$qname}){
                                    ${$target_info{$tname}}{$qname} = $align_info;
                                }
                            }
                        }
                    }
                    $pre_tend = $tend;
                    $pre_tstart = $tstart;
                    $pre_qend = $qend;
                    $pre_qstart = $qstart;
                    $pre_overlen = $overlen;
                    $pre_ident = $ident;
                }
                if ((($pre_qstart > 0) and ($pre_qstart <= 3) and ($pre_qend >= $qlen - 2)) or (($pre_tstart > 0) and ($pre_tstart <= 3) and ($pre_tend >= $tlen - 2))){
                    my $align_info = "$qname==$pre_qstart==$pre_qend==$tname==$pre_tstart==$pre_tend==$strand==$pre_overlen==$pre_ident";
                    if (!exists ${$target_info{$tname}}{$qname}){
                        ${$target_info{$tname}}{$qname} = $align_info;
                    }
                }
                elsif ((($pre_qstart > 0) and ($pre_qstart <= 3) and ($pre_tend >=  $tlen - 2) and ($strand eq 'Plus')) or (($pre_tstart > 0) and ($pre_tstart <= 3) and ($pre_qend >= $qlen - 2) and ($strand eq 'Plus')) or (($pre_qstart > 0) and ($pre_qstart <= 3) and ($pre_tstart <= 3) and ($strand eq 'Minus')) or (($pre_qend >= $qlen - 2) and ($pre_tend >=  $tlen - 2) and ($strand eq 'Minus'))){
                    my $align_info = "$qname==$pre_qstart==$pre_qend==$tname==$pre_tstart==$pre_tend==$strand==$pre_overlen==$pre_ident";
                    if (!exists ${$target_info{$tname}}{$qname}){
                        ${$target_info{$tname}}{$qname} = $align_info;
                    }
                }
            }
        }
    }
    %partial_match_5 = ();
    &filter_blast_result () if ($tag eq 'T');
    &filter_blast_result_Q () if ($tag eq 'Q');
    %target_info = ();
}

my $threshold_score = 0;
my $ave_log_score = 0;
my $ave_log_score_2 = 0;
my $standard_thresfold_score = 2.5;
$standard_thresfold_score += $adjust_score;

if ($PE_align_tag == 1){
    open (OUT, ">>$score_file");
    foreach my $align (keys %evaluate_pairs){
        my ($qname, $qstart, $qend, $tname, $tstart, $tend, $strand, $overlen, $identity) = split (/==/, $align);
        my $rate_PE = $evaluate_pairs{$align};
        my $score = 0;
        if ($rate_PE >= 0){
            $score = &validate_align ($overlen, $identity, $rate_PE);
            $evaluate_pairs{$align} = $score;
        }
        else{
            if (($overlen >= $min_match_length) and ($identity >= $min_identity)){
                $score = 1000;
                $evaluate_pairs{$align} = $score;
            }
            else{
                delete $evaluate_pairs{$align};
            }
        }
        print OUT "$qname\t$qstart\t$qend\t$tname\t$tstart\t$tend\t$strand\t$overlen\t$identity\t$rate_PE\t$score\n";
    }
    die "PE-reads appear not to be aligned\n" if ($num_score == 0);
    $ave_log_score = int (log ($sum_score / $num_score) * 100) / 100;
    $ave_log_score_2 = int (log ($sum_score_2 / $num_score_2) * 100) / 100 if ($num_score_2 > 0);
    if (($num_score_2 > 0) and ($ave_log_score + 0.5 < $ave_log_score_2)){
        $ave_log_score = $ave_log_score_2;
    }
    $threshold_score = $standard_thresfold_score - (5.5 - $ave_log_score);
    print OUT "##Ave log score: $ave_log_score Ave log score (rate_PE > 0): $ave_log_score_2\n";
    print OUT "##EOF";
    close (OUT);
}

if ($PE_align_tag == 1){
    print STDERR "\nThreshold score: $threshold_score (standard score: 2.5)\n\n";
    print OUTLOG "\nThreshold score: $threshold_score (standard score: 2.5)\n\n";
    foreach my $align (keys %evaluate_pairs){
        my $qunmatch_len = 0;
        my ($qname, $qstart, $qend, $tname, $tstart, $tend, $strand, $overlen, $identity) = split (/==/, $align);
        if (($evaluate_pairs{$align} >= $threshold_score) or (($overlen >= $min_match_length) and ($identity >= $min_identity))){
            my $Qlen = length $query_seq{$qname};
            my $Tlen = length $target_subcon_seq{$tname};
            $qunmatch_len = $Qlen - $qend if ($qstart <= 3) and ($qend + 2 < $Qlen);
            $qunmatch_len = $qstart - 1 if ($qstart > 3) and ($qend + 2 >= $Qlen);
            my $tname_base = $1 if ($tname =~/(\S+)\/(\d+)$/);
            my $tname_num = $2;
            $tname_num =~ s/^0*//;
            push @{${$TQ_cover{$tname_base}}{$qname}}, $tname_num;
            ${$TQ_strand{$tname_base}}{$qname} .= $strand . '=';
            ${$Q_match{$qname}}{$tname_base} += $qunmatch_len;
            if ($tstart <= 3){
                ${$T_match_5{$tname}}{$qname} = $align;
            }
            else{
                ${$T_match_3{$tname}}{$qname} = $align;
            }
        }
    }
}
else{
    foreach my $align (keys %evaluate_pairs){
        my $qunmatch_len = 0;
        my ($qname, $qstart, $qend, $tname, $tstart, $tend, $strand, $overlen, $identity) = split (/==/, $align);
        if (($overlen >= $min_match_length) and ($identity >= $min_identity)){
            my $Qlen = length $query_seq{$qname};
            my $Tlen = length $target_subcon_seq{$tname};
            $qunmatch_len = $Qlen - $qend if ($qstart <= 3) and ($qend + 2 < $Qlen);
            $qunmatch_len = $qstart - 1 if ($qstart > 3) and ($qend + 2 >= $Qlen);
            my $tname_base = $1 if ($tname =~/(\S+)\/(\d+)$/);
            my $tname_num = $2;
            $tname_num =~ s/^0*//;
            push @{${$TQ_cover{$tname_base}}{$qname}}, $tname_num;
            ${$TQ_strand{$tname_base}}{$qname} .= $strand . '=';
            ${$Q_match{$qname}}{$tname_base} += $qunmatch_len;
            if ($tstart <= 3){
                ${$T_match_5{$tname}}{$qname} = $align;
            }
            else{
                ${$T_match_3{$tname}}{$qname} = $align;
            }
        }
    }
}
%evaluate_pairs = ();

close (OUTLOG);

my %excluded_qname;

foreach my $qname (keys %Q_match){                  # remove target-query matches except for the best target-query pair with the longest overlapped sequences
    my $count = 0;
    foreach my $tname_base (sort {${$Q_match{$qname}}{$b} <=> ${$Q_match{$qname}}{$a}} keys %{$Q_match{$qname}}){
        $count++;
        if ($count > 1){
            delete ${$TQ_cover{$tname_base}}{$qname};
            delete ${$TQ_strand{$tname_base}}{$qname};
            my $match_flag = 0;
            foreach my $tname (sort keys %T_match_A){
                my $tname_b = $1 if ($tname =~/(\S+)\/(\d+)$/);
                my $num = $2;
                if ($tname_b eq $tname_base){
                    $match_flag = 1;
                    foreach my $qname2 (keys %{$T_match_A{$tname}}){
                        delete ${$T_match_A{$tname}}{$qname} if ($qname2 eq $qname);
                        $excluded_qname{$qname} = $tname_base;
                    }
                }
                elsif ($match_flag == 1){
                    last;
                }
            }
            $match_flag = 0;
            foreach my $tname (sort keys %T_match_5){
                my $tname_b = $1 if ($tname =~/(\S+)\/(\d+)$/);
                my $num = $2;
                if ($tname_b eq $tname_base){
                    $match_flag = 1;
                    foreach my $qname2 (keys %{$T_match_5{$tname}}){
                        delete ${$T_match_5{$tname}}{$qname} if ($qname2 eq $qname);
                        $excluded_qname{$qname} = $tname_base;
                    }
                }
                elsif ($match_flag == 1){
                    last;
                }
            }
            $match_flag = 0;
            foreach my $tname (sort keys %T_match_3){
                my $tname_b = $1 if ($tname =~/(\S+)\/(\d+)$/);
                my $num = $2;
                if ($tname_b eq $tname_base){
                    $match_flag = 1;
                    foreach my $qname2 (keys %{$T_match_3{$tname}}){
                        delete ${$T_match_3{$tname}}{$qname} if ($qname2 eq $qname);
                        $excluded_qname{$qname} = $tname_base;
                    }
                }
                elsif ($match_flag == 1){
                    last;
                }
            }
        }
    }
}
%Q_match = ();

my %skip_subcon;

foreach my $tname (keys %TQ_cover){             # remove query alignments with different alignment strand or discontinuous subseqeunces, out of the queries that cover multiple subsequences of a target scaffold
    foreach my $qname (keys %{$TQ_cover{$tname}}){
        my %count;
        @{${$TQ_cover{$tname}}{$qname}} = grep {!$count{$_}++} @{${$TQ_cover{$tname}}{$qname}};
        if (@{${$TQ_cover{$tname}}{$qname}} > 1){
            my $unallowed_subcon_flag = 0;       
            if ((${$TQ_strand{$tname}}{$qname} =~ /Plus/) and (${$TQ_strand{$tname}}{$qname} =~ /Minus/)){  # check wherher the strand direction of the aligned subseqeunces are the same (all Plus or all minus)
                $unallowed_subcon_flag = 1;
            }
            my $pre_num = 0;
            my $count = 0;
            foreach my $num (sort {$a <=> $b} @{${$TQ_cover{$tname}}{$qname}}){
                $count ++;
                if (($pre_num > 0) and ($unallowed_subcon_flag == 0)){         # check whether the aligned target subsequences are discontinuous
                    if ($num > $pre_num + 1){
                        for (my $i = $pre_num + 1; $i < $num; $i++){
                            my $subcon_num = sprintf ("%04d", $i);
                            my $subcon_name = $tname . '/' . $subcon_num;
                            if (length $target_subcon_seq{$subcon_name} > 200){
                                $unallowed_subcon_flag = 1;
                                last;
                            }
                            else{
                                $skip_subcon{$subcon_name} = $qname;
                            }
                        }
                    }
                }
                $pre_num = $num;
                if ($unallowed_subcon_flag == 0){                              # when there are >= 3 aligned subseqeunces, check whether all the region of the intervened subseqeunces are covered with a query
                    if ((@{${$TQ_cover{$tname}}{$qname}} >= 3) and ($count != 1) and ($count != @{${$TQ_cover{$tname}}{$qname}})){
                        my $subcon_num = sprintf ("%04d", $num);
                        my $subcon_name = $tname . '/' . $subcon_num;
                        if (!exists ${$T_match_A{$subcon_name}}{$qname}){
                            $unallowed_subcon_flag = 1;
                        }
                    }
                }
            }
            if ($unallowed_subcon_flag == 1){
                foreach my $num (@{${$TQ_cover{$tname}}{$qname}}){
                    my $subcon_num = sprintf ("%04d", $num);
                    my $subcon_name = $tname . '/' . $subcon_num;
                    delete ${$T_match_5{$subcon_name}}{$qname};
                    delete ${$T_match_3{$subcon_name}}{$qname};
                    delete ${$T_match_A{$subcon_name}}{$qname};
                }
                delete ${$TQ_cover{$tname}}{$qname};
                $excluded_qname{$qname} = $tname;
            }
        }
    }
}

foreach my $tname (keys %T_match_5){                                            # remove query alignments that have no overlap with other queries aligned to the same target terminus or remove query alignments with < $min_qalign_num
    my $qalign_num = scalar keys %{$T_match_5{$tname}};
    $qalign_num += scalar keys %{$T_match_A{$tname}} if (exists $T_match_A{$tname});
    if ($qalign_num >= 2){
        my %Qalign_num;
        foreach my $qname (keys %{$T_match_5{$tname}}){
            my $count_match = 1;
            foreach my $qname2 (keys %{$T_match_5{$tname}}){
                next if ($qname eq $qname2);
                if (exists ${$Q_match_5{$qname}}{$qname2} or exists ${$Q_match_5{$qname2}}{$qname} or exists ${$Q_match_3{$qname}}{$qname2} or ${$Q_match_3{$qname2}}{$qname} or exists ${$included_tname_Q{$qname}}{$qname2} or exists ${$included_tname_Q{$qname2}}{$qname} or exists ${$included_qname_Q{$qname}}{$qname2} or exists ${$included_qname_Q{$qname2}}{$qname}){
                    $count_match ++;
                }
            }
            if (exists $T_match_A{$tname}){
                foreach my $qname2 (keys %{$T_match_A{$tname}}){
                    next if ($qname eq $qname2);
                    if (exists ${$Q_match_5{$qname}}{$qname2} or exists ${$Q_match_5{$qname2}}{$qname} or exists ${$Q_match_3{$qname}}{$qname2} or ${$Q_match_3{$qname2}}{$qname} or exists ${$included_tname_Q{$qname}}{$qname2} or exists ${$included_tname_Q{$qname2}}{$qname} or exists ${$included_qname_Q{$qname}}{$qname2} or exists ${$included_qname_Q{$qname2}}{$qname}){
                        $count_match ++;
                    }
                }
            }
            $Qalign_num{$qname} = $count_match;
        }
        foreach my $qname3 (sort {$Qalign_num{$b} <=> $Qalign_num{$a}} keys %Qalign_num){
            my $num = $Qalign_num{$qname3};
            if (($num < $qalign_num / 2) or ($num < $min_qalign_num)){
                my $tname_base = $1 if ($tname =~ /(\S+)\/(\d+)$/);
                my $tname_num = $2;
                $tname_num =~ s/^0*//;
                delete ${$T_match_5{$tname}}{$qname3};
                if (@{${$TQ_cover{$tname_base}}{$qname3}} == 1){
                    delete ${$TQ_cover{$tname_base}}{$qname3};
                }
                else{
                    my $count = 0;
                    foreach (@{${$TQ_cover{$tname_base}}{$qname3}}){
                        last if ($_ == $tname_num);
                        $count ++;
                    }
                    splice (@{${$TQ_cover{$tname_base}}{$qname3}}, $count, 1) if ($count < @{${$TQ_cover{$tname_base}}{$qname3}});
                }
            }
        }
    }
    elsif ($qalign_num == 1){
        if ($min_qalign_num > 1){
            foreach my $qname (keys %{$T_match_5{$tname}}){
                my $tname_base = $1 if ($tname =~ /(\S+)\/(\d+)$/);
                my $tname_num = $2;
                $tname_num =~ s/^0*//;
                delete ${$T_match_5{$tname}}{$qname};
                if (@{${$TQ_cover{$tname_base}}{$qname}} == 1){
                    delete ${$TQ_cover{$tname_base}}{$qname};
                }
                else{
                    my $count = 0;
                    foreach (@{${$TQ_cover{$tname_base}}{$qname}}){
                        last if ($_ == $tname_num);
                        $count ++;
                    }
                    splice (@{${$TQ_cover{$tname_base}}{$qname}}, $count, 1) if ($count < @{${$TQ_cover{$tname_base}}{$qname}});
                }
            }
        }
    }
}

foreach my $tname (keys %T_match_3){                                            # remove query alignments that have no overlap with other queries aligned to the same target terminus or remove query alignments with < $min_qalign_num
    my $qalign_num = scalar keys %{$T_match_3{$tname}};
    $qalign_num += scalar keys %{$T_match_A{$tname}} if (exists $T_match_A{$tname});
    if ($qalign_num >= 2){
        my %Qalign_num;
        foreach my $qname (keys %{$T_match_3{$tname}}){
            my $count_match = 1;
            foreach my $qname2 (keys %{$T_match_3{$tname}}){
                next if ($qname eq $qname2);
                if (exists ${$Q_match_5{$qname}}{$qname2} or exists ${$Q_match_5{$qname2}}{$qname} or exists ${$Q_match_3{$qname}}{$qname2} or ${$Q_match_3{$qname2}}{$qname} or exists ${$included_tname_Q{$qname}}{$qname2} or exists ${$included_tname_Q{$qname2}}{$qname} or exists ${$included_qname_Q{$qname}}{$qname2} or exists ${$included_qname_Q{$qname2}}{$qname}){
                    $count_match ++;
                }
            }
            if (exists $T_match_A{$tname}){
                foreach my $qname2 (keys %{$T_match_A{$tname}}){
                    next if ($qname eq $qname2);
                    if (exists ${$Q_match_5{$qname}}{$qname2} or exists ${$Q_match_5{$qname2}}{$qname} or exists ${$Q_match_3{$qname}}{$qname2} or ${$Q_match_3{$qname2}}{$qname} or exists ${$included_tname_Q{$qname}}{$qname2} or exists ${$included_tname_Q{$qname2}}{$qname} or exists ${$included_qname_Q{$qname}}{$qname2} or exists ${$included_qname_Q{$qname2}}{$qname}){
                        $count_match ++;
                    }
                }
            }
            $Qalign_num{$qname} = $count_match;
        }
        foreach my $qname3 (sort {$Qalign_num{$b} <=> $Qalign_num{$a}} keys %Qalign_num){
            my $num = $Qalign_num{$qname3};
            if (($num < $qalign_num / 2) or ($num < $min_qalign_num)){
                my $tname_base = $1 if ($tname =~ /(\S+)\/(\d+)$/);
                my $tname_num = $2;
                $tname_num =~ s/^0*//;
                delete ${$T_match_3{$tname}}{$qname3};
                if (@{${$TQ_cover{$tname_base}}{$qname3}} == 1){
                    delete ${$TQ_cover{$tname_base}}{$qname3};
                }
                else{
                    my $count = 0;
                    foreach (@{${$TQ_cover{$tname_base}}{$qname3}}){
                        last if ($_ == $tname_num);
                        $count ++;
                    }
                    splice (@{${$TQ_cover{$tname_base}}{$qname3}}, $count, 1) if ($count < @{${$TQ_cover{$tname_base}}{$qname3}});
                }
            }
        }
    }
    elsif ($qalign_num == 1){
        if ($min_qalign_num > 1){
            foreach my $qname (keys %{$T_match_3{$tname}}){
                my $tname_base = $1 if ($tname =~ /(\S+)\/(\d+)$/);
                my $tname_num = $2;
                $tname_num =~ s/^0*//;
                delete ${$T_match_3{$tname}}{$qname};
                if (@{${$TQ_cover{$tname_base}}{$qname}} == 1){
                    delete ${$TQ_cover{$tname_base}}{$qname};
                }
                else{
                    my $count = 0;
                    foreach (@{${$TQ_cover{$tname_base}}{$qname}}){
                        last if ($_ == $tname_num);
                        $count ++;
                    }
                    splice (@{${$TQ_cover{$tname_base}}{$qname}}, $count, 1) if ($count < @{${$TQ_cover{$tname_base}}{$qname}});
                }
            }
        }
    }
}
%Q_match_5 = ();
%Q_match_3 = ();
%included_tname_Q = ();

foreach my $tname (sort keys %T_match_5){
    foreach my $qname (sort keys %{$T_match_5{$tname}}){
        ${$T_match_ALL{$tname}}{$qname} = ${$T_match_5{$tname}}{$qname};
    }
}

foreach my $tname (sort keys %T_match_3){
    foreach my $qname (sort keys %{$T_match_3{$tname}}){
        ${$T_match_ALL{$tname}}{$qname} = ${$T_match_3{$tname}}{$qname};
    }
}

foreach my $tname (sort keys %T_match_A){
    foreach my $qname (sort keys %{$T_match_A{$tname}}){
        ${$T_match_ALL{$tname}}{$qname} = ${$T_match_A{$tname}}{$qname};
    }
}

foreach my $tname (sort keys %T_match_ALL){         # filterng overlapped alignments with only a single overlapping subcontig
    my $num_query = scalar keys %{$T_match_ALL{$tname}};
    my %multi_query;
    my %count_5;
    my %count_3;
    my %count_A;
    my $base_name = $1 if ($tname =~/(\S+)\/(\d+)$/);
    my $num = $2;
    if ($num_query > 1){
        foreach my $qname (sort keys %{$T_match_ALL{$tname}}){
            if (exists ${$T_match_A{$tname}}{$qname}){
                $count_A{$qname} = 1;
            }
            elsif (exists ${$T_match_5{$tname}}{$qname}){
                $count_5{$qname} = 1;
            }
            elsif (exists ${$T_match_3{$tname}}{$qname}){
                $count_3{$qname} = 1;
            }
            $multi_query{$qname} = 1;
        }
        if (scalar keys %count_A > 0){
            foreach my $mqname (%multi_query){
                if (!exists $count_A{$mqname}){
                    if ((exists ${$TQ_cover{$base_name}}{$mqname}) and (@{${$TQ_cover{$base_name}}{$mqname}} == 1)){
                        delete ${$T_match_ALL{$tname}}{$mqname};
                        delete ${$TQ_cover{$base_name}}{$mqname};
                        $excluded_qname{$mqname} = $tname;
                    }
                }
            }
        }
        elsif (scalar keys %count_5 > 1){
            my %match_single;
            my %match_multi;
            foreach my $qname5 (%count_5){
                next if (!exists ${$TQ_cover{$base_name}}{$qname5});
                my $match_num = scalar @{${$TQ_cover{$base_name}}{$qname5}};
                $match_single{$qname5} = 1 if ($match_num == 1);
                $match_multi{$qname5} = 1 if ($match_num > 1);
            }
            if ((scalar keys %match_single > 0) and (scalar keys %match_multi > 0)){
                foreach my $qname5 (%match_single){
                    delete ${$T_match_ALL{$tname}}{$qname5};
                    delete ${$TQ_cover{$base_name}}{$qname5};
                    $excluded_qname{$qname5} = $tname;
                }
            }
            elsif ((scalar keys %match_single > 1) and (scalar keys %match_multi == 0)){
                foreach my $qname5 (%match_single){
                    if (exists ${$T_match_ALL{$tname}}{$qname5}){
                        my ($qname1, $qstart, $qend, $tname1, $tstart, $tend, $strand, $overlen) = split (/==/, ${$T_match_ALL{$tname}}{$qname5});
                        my $unmatch_len = $qstart - 1 + length ($query_seq{$qname1}) - $qend;
                        my $highest_match_query = '';
                        my $highest_unmatch_len = 0;
                        if ($highest_unmatch_len == 0){
                            $highest_match_query = $qname5;
                            $highest_unmatch_len = $unmatch_len;
                            next;
                        }
                        if ($highest_unmatch_len > $unmatch_len){
                            delete ${$T_match_ALL{$tname}}{$qname5};
                            delete ${$TQ_cover{$base_name}}{$qname5};
                            $excluded_qname{$qname5} = $tname;
                        }
                        else{
                            delete ${$T_match_ALL{$tname}}{$highest_match_query};
                            delete ${$TQ_cover{$base_name}}{$highest_match_query};
                            $excluded_qname{$highest_match_query} = $tname;
                            $highest_match_query = $qname5;
                            $highest_unmatch_len = $unmatch_len;
                        }
                    }
                }
            }
        }
        elsif (scalar keys %count_3 > 1){
            my %match_single;
            my %match_multi;
            foreach my $qname3 (%count_3){
                next if (!exists ${$TQ_cover{$base_name}}{$qname3});
                my $match_num = scalar @{${$TQ_cover{$base_name}}{$qname3}};
                $match_single{$qname3} = 1 if ($match_num == 1);
                $match_multi{$qname3} = 1 if ($match_num > 1);
            }
            if ((scalar keys %match_single > 0) and (scalar keys %match_multi > 0)){
                foreach my $qname3 (%match_single){
                    delete ${$T_match_ALL{$tname}}{$qname3};
                    delete ${$TQ_cover{$base_name}}{$qname3};
                    $excluded_qname{$qname3} = $tname;
                }
            }
            elsif ((scalar keys %match_single > 1) and (scalar keys %match_multi == 0)){
                foreach my $qname3 (%match_single){
                    if (exists ${$T_match_ALL{$tname}}{$qname3}){
                        my ($qname1, $qstart, $qend, $tname1, $tstart, $tend, $strand, $overlen) = split (/==/, ${$T_match_ALL{$tname}}{$qname3});
                        my $unmatch_len = $qstart - 1 + length ($query_seq{$qname1}) - $qend;
                        my $highest_match_query = '';
                        my $highest_unmatch_len = 0;
                        if ($highest_unmatch_len == 0){
                            $highest_match_query = $qname3;
                            $highest_unmatch_len = $unmatch_len;
                            next;
                        }
                        if ($highest_unmatch_len > $unmatch_len){
                            delete ${$T_match_ALL{$tname}}{$qname3};
                            delete ${$TQ_cover{$base_name}}{$qname3};
                            $excluded_qname{$qname3} = $tname;
                        }
                        else{
                            delete ${$T_match_ALL{$tname}}{$highest_match_query};
                            delete ${$TQ_cover{$base_name}}{$highest_match_query};
                            $excluded_qname{$highest_match_query} = $tname;
                            $highest_match_query = $qname3;
                            $highest_unmatch_len = $unmatch_len;
                        }
                    }
                }
            }
        }
    }
}

my %gap_5;
my %gap_3;
my %gap_A;
my %term5;
my %term3;
my %used_longread_T5;
my %used_longread_T3;
my %used_longread_TA;

my $out_align_list = "$ite_dir$out_prefix-3m.selected-align-list.txt";
open (OUT1, "> $out_align_list");
print OUT1 '-' x 90, "\n";
print OUT1 "query_name\tpos1\tpos2\ttarget_name\tpos1\tpos2\tstrand\toverlap_len\tidentity\n";
print OUT1 '-' x 90, "\n";
foreach my $tname (sort keys %T_match_ALL){         # set query coordinates to gap regions
    next if (scalar keys %{$T_match_ALL{$tname}} == 0);
    my $base_name = $1 if ($tname =~/(\S+)\/(\d+)$/);
    my $num = $2;
    $num =~ s/^0*//;
    my $up_tname = $base_name . '/' . (sprintf ("%04d", $num - 1));
    my $down_tname = $base_name . '/' . (sprintf ("%04d", $num + 1));
    my $hit_flag = 0;
    my $up_tname_2 = $up_tname;
    my $down_tname_2 = $down_tname;
    if (($num >= 3) and (exists $skip_subcon{$up_tname}) and (exists $T_match_5{$tname} or exists $T_match_A{$tname})){
        while ($hit_flag == 0){
            $num --;
            $up_tname_2 = $base_name . '/' . (sprintf ("%04d", $num - 1));
            $hit_flag = 1 if (!exists $skip_subcon{$up_tname_2} and (exists $T_match_3{$up_tname_2} or exists $T_match_A{$up_tname_2}));
            last if ($num == 1);
        }
    }
    $up_tname = $up_tname_2 if ($hit_flag == 1);
    $hit_flag = 0;
    if (($num <= @{$Tsubcon_info{$base_name}} - 2) and (exists $skip_subcon{$down_tname}) and (exists $T_match_3{$tname} or exists $T_match_A{$tname})){
        while ($hit_flag == 0){
            $num ++;
            $down_tname_2 = $base_name . '/' . (sprintf ("%04d", $num + 1));
            $hit_flag = 1 if (!exists $skip_subcon{$down_tname_2} and (exists $T_match_5{$down_tname_2} or exists $T_match_A{$down_tname_2}));
            last if ($num == @{$Tsubcon_info{$base_name}});
        }
    }
    $down_tname = $down_tname_2 if ($hit_flag == 1);
    
    foreach my $qname (sort keys %{$T_match_ALL{$tname}}){
        my @out = split (/==/, ${$T_match_ALL{$tname}}{$qname});
        print OUT1 join ("\t", @out), "\n";
        my $Qlen = length $query_seq{$qname};
        my ($qname1, $qstart, $qend, $tname1, $tstart, $tend, $strand, $overlen) = split (/==/, ${$T_match_ALL{$tname}}{$qname});
        if (($extend == 1) and ($num == 1) and ($ite_num == 1)){
            if ((exists ${$T_match_5{$tname}}{$qname}) or (exists ${$T_match_A{$tname}}{$qname})){
                my $ex_start = 1 if ($strand eq 'Plus');
                $ex_start = $qend + 1 if ($strand eq 'Minus');
                my $ex_end = $qstart - 1 if ($strand eq 'Plus');
                $ex_end = $Qlen if ($strand eq 'Minus');
                my $ex_len = $ex_end - $ex_start + 1;
                if ($ex_len > 0){
                    if (exists $term5{$base_name}){
                        my ($qname2, $qstart2, $qend2, $extend_len2, $strand2) = split (/==/, $term5{$base_name});
                        if ($extend_len2 < $ex_len){
                            $term5{$base_name} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                            $used_longread_T5{$tname} = $qname;
                        }
                    }
                    else{
                        $term5{$base_name} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                        $used_longread_T5{$tname} = $qname;
                    }
                }
            }
        }
        elsif (($extend == 1) and ($num == @{$Tsubcon_info{$base_name}}) and ($ite_num == 1)){
            if ((exists ${$T_match_3{$tname}}{$qname}) or (exists ${$T_match_A{$tname}}{$qname})){
                my $ex_start = $qend + 1 if ($strand eq 'Plus');
                $ex_start = 1 if ($strand eq 'Minus');
                my $ex_end = $Qlen if ($strand eq 'Plus');
                $ex_end = $qstart - 1 if ($strand eq 'Minus');
                my $ex_len = $ex_end - $ex_start + 1;
                if ($ex_len > 0){
                    if (exists $term3{$base_name}){
                        my ($qname2, $qstart2, $qend2, $extend_len2, $strand2) = split (/==/, $term3{$base_name});
                        if ($extend_len2 < $ex_len){
                            $term3{$base_name} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                            $used_longread_T3{$tname} = $qname;
                        }
                    }
                    else{
                        $term3{$base_name} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                        $used_longread_T3{$tname} = $qname;
                    }
                }
            }
        }
        if (exists ${$T_match_A{$tname}}{$qname}){
            if ((!exists $gap_A{$up_tname}) and ($num != 1)){
                if (exists ${$T_match_ALL{$up_tname}}{$qname} and ((exists ${$T_match_A{$up_tname}}{$qname}) or (exists ${$T_match_3{$up_tname}}{$qname}))){
                    my ($u_qname1, $u_qstart, $u_qend, $u_tname1, $u_tstart, $u_tend, $u_strand, $u_overlen) = split (/==/, ${$T_match_ALL{$up_tname}}{$qname});
                    my $ex_start = $u_qend + 1 if ($strand eq 'Plus');
                    $ex_start = $qend + 1 if ($strand eq 'Minus');
                    my $ex_end = $qstart - 1 if ($strand eq 'Plus');
                    $ex_end = $u_qstart - 1 if ($strand eq 'Minus');
                    my $ex_len = $ex_end - $ex_start + 1;
                    $gap_A{$up_tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                }
                else{
                    my $ex_start = 1 if ($strand eq 'Plus');
                    $ex_start = $qend + 1 if ($strand eq 'Minus');
                    my $ex_end = $qstart - 1 if ($strand eq 'Plus');
                    $ex_end = $Qlen if ($strand eq 'Minus');
                    my $ex_len = $ex_end - $ex_start + 1;
                    next if ($ex_len < 1);
                    if (exists $gap_3{$up_tname}){
                        my ($qname2, $qstart2, $qend2, $extend_len2, $strand2) = split (/==/, $gap_3{$up_tname});
                        if ($extend_len2 < $ex_len){
                            $gap_3{$up_tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                        }
                    }
                    else{
                        $gap_3{$up_tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                    }
                }
            }
            if ((!exists $gap_A{$tname}) and ($num != @{$Tsubcon_info{$base_name}})){
                if (exists ${$T_match_ALL{$down_tname}}{$qname} and ((exists ${$T_match_A{$down_tname}}{$qname}) or (exists ${$T_match_5{$down_tname}}{$qname}))){
                    my ($d_qname1, $d_qstart, $d_qend, $d_tname1, $d_tstart, $d_tend, $d_strand, $d_overlen) = split (/==/, ${$T_match_ALL{$down_tname}}{$qname});
                    my $ex_start = $qend + 1 if ($strand eq 'Plus');
                    $ex_start = $d_qend + 1 if ($strand eq 'Minus');
                    my $ex_end = $d_qstart - 1 if ($strand eq 'Plus');
                    $ex_end = $qstart - 1 if ($strand eq 'Minus');
                    my $ex_len = $ex_end - $ex_start + 1;
                    $gap_A{$tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                }
                else{
                    my $ex_start = $qend + 1 if ($strand eq 'Plus');
                    $ex_start = 1 if ($strand eq 'Minus');
                    my $ex_end = $Qlen if ($strand eq 'Plus');
                    $ex_end = $qstart - 1 if ($strand eq 'Minus');
                    my $ex_len = $ex_end - $ex_start + 1;
                    next if ($ex_len < 1);
                    if (exists $gap_5{$tname}){
                        my ($qname2, $qstart2, $qend2, $extend_len2, $strand2) = split (/==/, $gap_5{$tname});
                        if ($extend_len2 < $ex_len){
                            $gap_5{$tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                        }
                    }
                    else{
                        $gap_5{$tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                    }
                }
            }
        }
        elsif ((exists ${$T_match_5{$tname}}{$qname}) and (!exists $gap_A{$up_tname}) and ($num != 1)){
            if (exists ${$T_match_ALL{$up_tname}}{$qname} and exists ${$T_match_3{$up_tname}}{$qname}){
                my ($u_qname1, $u_qstart, $u_qend, $u_tname1, $u_tstart, $u_tend, $u_strand, $u_overlen) = split (/==/, ${$T_match_ALL{$up_tname}}{$qname});
                my $ex_start = $u_qend + 1 if ($strand eq 'Plus');
                $ex_start = $qend + 1 if ($strand eq 'Minus');
                my $ex_end = $qstart - 1 if ($strand eq 'Plus');
                $ex_end = $u_qstart - 1 if ($strand eq 'Minus');
                my $ex_len = $ex_end - $ex_start + 1;
                $gap_A{$up_tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
            }
            else{
                my $ex_start = 1 if ($strand eq 'Plus');
                $ex_start = $qend + 1 if ($strand eq 'Minus');
                my $ex_end = $qstart - 1 if ($strand eq 'Plus');
                $ex_end = $Qlen if ($strand eq 'Minus');
                my $ex_len = $ex_end - $ex_start + 1;
                next if ($ex_len < 1);
                if (exists $gap_3{$up_tname}){
                    my ($qname2, $qstart2, $qend2, $extend_len2, $strand2) = split (/==/, $gap_3{$up_tname});
                    if ($extend_len2 < $ex_len){
                        $gap_3{$up_tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                    }
                }
                else{
                    $gap_3{$up_tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                }
            }
        }
        elsif ((exists ${$T_match_3{$tname}}{$qname}) and (!exists $gap_A{$tname}) and ($num != @{$Tsubcon_info{$base_name}})){
            if (exists ${$T_match_ALL{$down_tname}}{$qname} and exists ${$T_match_5{$down_tname}}{$qname}){
                my ($d_qname1, $d_qstart, $d_qend, $d_tname1, $d_tstart, $d_tend, $d_strand, $d_overlen) = split (/==/, ${$T_match_ALL{$down_tname}}{$qname});
                my $ex_start = $qend + 1 if ($strand eq 'Plus');
                $ex_start = $d_qend + 1 if ($strand eq 'Minus');
                my $ex_end = $d_qstart - 1 if ($strand eq 'Plus');
                $ex_end = $qstart - 1 if ($strand eq 'Minus');
                my $ex_len = $ex_end - $ex_start + 1;
                $gap_A{$tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
            }
            else{
                my $ex_start = $qend + 1 if ($strand eq 'Plus');
                $ex_start = 1 if ($strand eq 'Minus');
                my $ex_end = $Qlen if ($strand eq 'Plus');
                $ex_end = $qstart - 1 if ($strand eq 'Minus');
                my $ex_len = $ex_end - $ex_start + 1;
                next if ($ex_len < 1);
                if (exists $gap_5{$tname}){
                    my ($qname2, $qstart2, $qend2, $extend_len2, $strand2) = split (/==/, $gap_5{$tname});
                    if ($extend_len2 < $ex_len){
                        $gap_5{$tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                    }
                }
                else{
                    $gap_5{$tname} = $qname . '==' . $ex_start . '==' . $ex_end . '==' . $ex_len . '==' . $strand;
                }
            }
        }
    }
}
close (OUT1);

undef %T_match_5;
undef %T_match_3;
undef %T_match_A;
undef %T_match_I;
undef %T_match_ALL;
undef %TQ_cover;
undef %TQ_strand;
undef %included_tname;
undef %included_qname_Q;


my %gap_Aseq;
my %gap_Cseq;
my @subcon_subfiles;
my @T_PEalignF_subfiles;
my @T_PEalignR_subfiles;
my @gapA_subfiles;
my @gap5_subfiles;
my @gap3_subfiles;
my @connect_subfiles;
my @delete_subfiles;
if (($thread_connect > 1) and ($Config{useithreads})){          # connect neighboring subcontigs
    $localtime = localtime;
    print STDERR "\nContig connection start in parallel: $localtime\n";
    my $scaf_num = scalar keys %target_seq;
    my $no_scaf_num = scalar keys %target_noscaf;
    my $seq_num = $scaf_num - $no_scaf_num;
    $thread_connect = 24 if ($thread_connect > 24);
    my $div_num = int ($seq_num / $thread_connect) + 1;
    my $count_sub = 1;
    my $count_tbase = 0;
    my @target_subfiles;
    my $pre_tbase = '';
    my %tsubseq;
    foreach my $tname (sort keys %target_subcon_seq){
        my $tname_base = $1 if ($tname =~/(\S+)\/(\d+)$/);
        my $tname_num = $2;
        next if (exists $target_noscaf{$tname_base});
        $pre_tbase = $tname_base if ($pre_tbase eq '');
        if ($pre_tbase ne $tname_base){
            $count_tbase ++;
            if ($count_tbase >= $div_num){
                my $target_subcon_seq_tempout = "$temp_dir/$out_prefix-subcon-seq-$count_sub.fa";
                my $target_PEalignF_tempout = "$temp_dir/$out_prefix-target-PE-align-F-$count_sub.txt";
                my $target_PEalignR_tempout = "$temp_dir/$out_prefix-target-PE-align-R-$count_sub.txt";
                my $target_gapA_tempout = "$temp_dir/$out_prefix-gapA-$count_sub.txt";
                my $target_gap5_tempout = "$temp_dir/$out_prefix-gap5-$count_sub.txt";
                my $target_gap3_tempout = "$temp_dir/$out_prefix-gap3-$count_sub.txt";
                open (OUT, "> $target_subcon_seq_tempout");
                open (OUT2, "> $target_PEalignF_tempout") if ($PE_align_tag == 1);
                open (OUT3, "> $target_PEalignR_tempout") if ($PE_align_tag == 1);
                open (OUT4, "> $target_gapA_tempout");
                open (OUT5, "> $target_gap5_tempout");
                open (OUT6, "> $target_gap3_tempout");
                foreach my $tname_2 (sort keys %tsubseq){
                    print OUT ">$tname_2 ";
                    if (exists $target_gap_size{$tname_2}){
                        print OUT "$target_gap_size{$tname_2}\n";
                    }
                    else{
                        print OUT "\n";
                    }
                    print OUT $tsubseq{$tname_2}, "\n";
                    print OUT4 "$tname_2\t", $gap_A{$tname_2}, "\n" if (exists $gap_A{$tname_2});
                    print OUT5 "$tname_2\t", $gap_5{$tname_2}, "\n" if (exists $gap_5{$tname_2});
                    print OUT6 "$tname_2\t", $gap_3{$tname_2}, "\n" if (exists $gap_3{$tname_2});
                    if ($PE_align_tag == 1){
                        foreach my $pos (keys %{$target_read_pos_F{$tname_2}}){
                            print OUT2 "$tname_2\t$pos\t${$target_read_pos_F{$tname_2}}{$pos}\n";
                        }
                        foreach my $pos (keys %{$target_read_pos_R{$tname_2}}){
                            print OUT3 "$tname_2\t$pos\t${$target_read_pos_R{$tname_2}}{$pos}\n";
                        }
                    }
                }
                close (OUT);
                close (OUT2) if ($PE_align_tag == 1);
                close (OUT3) if ($PE_align_tag == 1);
                close (OUT4);
                close (OUT5);
                close (OUT6);
                push @subcon_subfiles, $target_subcon_seq_tempout;
                push @T_PEalignF_subfiles, $target_PEalignF_tempout if ($PE_align_tag == 1);
                push @T_PEalignR_subfiles, $target_PEalignR_tempout if ($PE_align_tag == 1);
                push @gapA_subfiles, $target_gapA_tempout;
                push @gap5_subfiles, $target_gap5_tempout;
                push @gap3_subfiles, $target_gap3_tempout;
                $count_sub ++;
                $count_tbase = 0;
                %tsubseq = ();
            }
        }
        $tsubseq{$tname} = $target_subcon_seq{$tname};
        $pre_tbase = $tname_base;
    }
    my $query_PEalignF_tempout = "$temp_dir/$out_prefix-query-PE-align-F.txt";
    my $query_PEalignR_tempout = "$temp_dir/$out_prefix-query-PE-align-R.txt";
    if ($PE_align_tag == 1){
        open (OUT2, "> $query_PEalignF_tempout") if ($PE_align_tag == 1);
        open (OUT3, "> $query_PEalignR_tempout") if ($PE_align_tag == 1);
        foreach my $qname (keys %query_read_pos_F){
            foreach my $pos (keys %{$query_read_pos_F{$qname}}){
                print OUT2 "$qname\t$pos\t${$query_read_pos_F{$qname}}{$pos}\n";
            }
            foreach my $pos (keys %{$query_read_pos_R{$qname}}){
                print OUT3 "$qname\t$pos\t${$query_read_pos_R{$qname}}{$pos}\n";
            }
        }
        close (OUT2) if ($PE_align_tag == 1);
        close (OUT3) if ($PE_align_tag == 1);  
    }
    undef %query_read_pos_F;
    undef %query_read_pos_R;
    
    $count_tbase ++;
    my $target_subcon_seq_tempout = "$temp_dir/$out_prefix-subcon-seq-$count_sub.fa";
    my $target_PEalignF_tempout = "$temp_dir/$out_prefix-PE-align-F-$count_sub.txt";
    my $target_PEalignR_tempout = "$temp_dir/$out_prefix-PE-align-R-$count_sub.txt";
    my $target_gapA_tempout = "$temp_dir/$out_prefix-gapA-$count_sub.txt";
    my $target_gap5_tempout = "$temp_dir/$out_prefix-gap5-$count_sub.txt";
    my $target_gap3_tempout = "$temp_dir/$out_prefix-gap3-$count_sub.txt";
    open (OUT, "> $target_subcon_seq_tempout");
    open (OUT2, "> $target_PEalignF_tempout") if ($PE_align_tag == 1);
    open (OUT3, "> $target_PEalignR_tempout") if ($PE_align_tag == 1);
    open (OUT4, "> $target_gapA_tempout");
    open (OUT5, "> $target_gap5_tempout");
    open (OUT6, "> $target_gap3_tempout");
    foreach my $tname_2 (sort keys %tsubseq){
        print OUT ">$tname_2 ";
        if (exists $target_gap_size{$tname_2}){
            print OUT "$target_gap_size{$tname_2}\n";
        }
        else{
            print OUT "\n";
        }
        print OUT $tsubseq{$tname_2}, "\n";
        print OUT4 "$tname_2\t", $gap_A{$tname_2}, "\n" if (exists $gap_A{$tname_2});
        print OUT5 "$tname_2\t", $gap_5{$tname_2}, "\n" if (exists $gap_5{$tname_2});
        print OUT6 "$tname_2\t", $gap_3{$tname_2}, "\n" if (exists $gap_3{$tname_2});
        if ($PE_align_tag == 1){
            foreach my $pos (keys %{$target_read_pos_F{$tname_2}}){
                print OUT2 "$tname_2\t$pos\t${$target_read_pos_F{$tname_2}}{$pos}\n";
            }
            foreach my $pos (keys %{$target_read_pos_R{$tname_2}}){
                print OUT3 "$tname_2\t$pos\t${$target_read_pos_R{$tname_2}}{$pos}\n";
            }
        }
    }
    close (OUT);
    close (OUT2) if ($PE_align_tag == 1);
    close (OUT3) if ($PE_align_tag == 1);
    close (OUT4);
    close (OUT5);
    close (OUT6);
    push @subcon_subfiles, $target_subcon_seq_tempout;
    push @T_PEalignF_subfiles, $target_PEalignF_tempout if ($PE_align_tag == 1);
    push @T_PEalignR_subfiles, $target_PEalignR_tempout if ($PE_align_tag == 1);
    push @gapA_subfiles, $target_gapA_tempout;
    push @gap5_subfiles, $target_gap5_tempout;
    push @gap3_subfiles, $target_gap3_tempout;
    undef %tsubseq;
    undef %target_read_pos_F;
    undef %target_read_pos_R;
    
    my @connect_jobs;
    my $num = 0;
    foreach my $subfile (@subcon_subfiles){
        $num ++;
        my $T_PEalignF_subfile = $T_PEalignF_subfiles[$num - 1] if ($PE_align_tag == 1);
        my $T_PEalignR_subfile = $T_PEalignR_subfiles[$num - 1] if ($PE_align_tag == 1);
        my $gapA_subfile = $gapA_subfiles[$num - 1];
        my $gap5_subfile = $gap5_subfiles[$num - 1];
        my $gap3_subfile = $gap3_subfiles[$num - 1];
        my ($thread) = threads->new(\&connect_subcon, $subfile, $num, $out_prefix, $temp_dir, $query_file, $gapA_subfile, $gap5_subfile, $gap3_subfile, $PE_align_tag, $max_insert_size, $read_length, $T_PEalignF_subfile, $T_PEalignR_subfile, $query_PEalignF_tempout, $query_PEalignR_tempout) if ($PE_align_tag == 1);
        ($thread) = threads->new(\&connect_subcon, $subfile, $num, $out_prefix, $temp_dir, $query_file, $gapA_subfile, $gap5_subfile, $gap3_subfile, $PE_align_tag) if ($PE_align_tag == 0);
        push @connect_jobs, $thread;
    }
    foreach (@connect_jobs){
        my ($connect_file, $delete_tqname_file, $subfile_num) = $_->join;
        print STDERR "Contig connection completed ($subfile_num)\n";
        push @connect_subfiles, $connect_file;
        push @delete_subfiles, $delete_tqname_file;
    }
    if (@connect_subfiles > 0){
        foreach my $subfile (@connect_subfiles){
            unless (-e $subfile){
                print STDERR "Caution: $subfile was not generated\n";
                next;
            }
            open (FILE, $subfile);
            while (my $line = <FILE>){
                chomp $line;
                next if ($line =~ /^$/);
                my ($tname, $tag, $gap_seq, $qname_set) = split (/\t/, $line);
                if ($tag eq 'A'){
                    $gap_Aseq{$tname} = $gap_seq;
                    $used_longread_TA{$tname} = $qname_set;
                }
                elsif ($tag eq 'C'){
                    $gap_Cseq{$tname} = $gap_seq;
                }
            }
            close (FILE);
        }
        foreach my $subfile (@delete_subfiles){
            unless (-e $subfile){
                print STDERR "Caution: $subfile was not generated\n";
                next;
            }
            open (FILE, $subfile);
            while (my $line = <FILE>){
                chomp $line;
                next if ($line =~ /^$/);
                my ($tname, $qname) = split (/\t/, $line);
                delete $gap_5{$tname};
                delete $gap_3{$tname};
            }
            close (FILE);
        }
    }
    undef @subcon_subfiles;
    undef @connect_jobs;
    undef @connect_subfiles;
    undef @delete_subfiles;
    system ("rm -f $temp_dir/$out_prefix-subcon-seq-[0-9].fa");
    system ("rm -f $temp_dir/$out_prefix-target-PE-align-[FR]-[0-9].txt");
    system ("rm -f $temp_dir/$out_prefix-query-PE-align-[FR].txt");
    system ("rm -f $temp_dir/$out_prefix-[0-9]-seq[12].fa");
    system ("rm -f $temp_dir/$out_prefix-gap[A53]-[0-9].txt");
    my @file_list = <$temp_dir/$out_prefix-[0-9][0-9]-seq[12].fa>;
    system ("rm -f $temp_dir/$out_prefix-[0-9][0-9]-seq[12].fa") if (@file_list > 0);
    system ("rm -f $temp_dir/$out_prefix-subcon-connect-[0-9]*");
    system ("rm -f $temp_dir/$out_prefix-delete_tqname-[0-9]*");
}
else{
    $localtime = localtime;
    print STDERR "\nContig connection start with single-thread: $localtime\n";
    foreach my $tname (sort keys %target_subcon_seq){               # merge overlaps of gap 5'- and 3'-terminal seqeunces of subcontigs surrounding filled gaps to complete a gap seqeunce
        my $base_name = $1 if ($tname =~/(\S+)\/(\d+)$/);
        my $num = $2;
        next if (exists $target_noscaf{$base_name});
        $num =~ s/^0*//;
        next if ($num == @{$Tsubcon_info{$base_name}});
        my $down_tname = $base_name . '/' . (sprintf ("%04d", $num + 1));
        my $gap_5_unmatch_len = 0;
        my $gap_3_unmatch_len = 0;
        my $gap5_seq = '';
        my $gap3_seq = '';
        my $gap5_qname = '';
        my $gap3_qname = '';
        my $gap5_TQtag = '';
        my $gap3_TQtag = '';
        if (!exists $gap_A{$tname}){
            if ((exists $gap_5{$tname}) or (exists $gap_3{$tname})){
                if (exists $gap_5{$tname}){
                    my ($qname5, $start5, $end5, $exlen5, $strand5) = split (/==/, $gap_5{$tname});
                    $gap_5_unmatch_len = $exlen5;
                    $gap5_seq = $query_seq{$qname5};
                    $gap5_seq = reverse $gap5_seq if ($strand5 eq 'Minus');
                    $gap5_seq =~ tr/ACGT/TGCA/ if ($strand5 eq 'Minus');
                    $gap5_qname = $qname5;
                    $gap5_TQtag = 'Q3' if ($strand5 eq 'Plus');
                    $gap5_TQtag = 'Q5' if ($strand5 eq 'Minus');
                }
                else{
                    $gap5_seq = $target_subcon_seq{$tname};
                    $gap5_TQtag = 'T3';
                }
                if (exists $gap_3{$tname}){
                    my ($qname3, $start3, $end3, $exlen3, $strand3) = split (/==/, $gap_3{$tname});
                    $gap_3_unmatch_len = $exlen3;
                    $gap3_seq = $query_seq{$qname3};
                    $gap3_seq = reverse $gap3_seq if ($strand3 eq 'Minus');
                    $gap3_seq =~ tr/ACGT/TGCA/ if ($strand3 eq 'Minus');
                    $gap3_qname = $qname3;
                    $gap3_TQtag = 'Q5' if ($strand3 eq 'Plus');
                    $gap3_TQtag = 'Q3' if ($strand3 eq 'Minus');
                }
                else{
                    $gap3_seq = $target_subcon_seq{$down_tname};
                    $gap3_TQtag = 'T5';
                }
                my ($match, $match_len1, $match_len2, $identity2) = &local_align3 (\$gap5_seq, \$gap3_seq);
                my $rate_PE = 0;
                my $read_num = 0;
                if (($match == 1) and ($PE_align_tag == 1)){
                    my $gap3_len = length $gap3_seq;
                    my $gap5_len = length $gap5_seq;
                    next if ($gap5_len < $read_length) or ($gap3_len < $read_length);
                    my $gap5_name1 = $tname if ($gap5_TQtag eq 'T3');
                    $gap5_name1 = $gap5_qname if ($gap5_TQtag eq 'Q5') or ($gap5_TQtag eq 'Q3');
                    my $gap3_name1 = $down_tname if ($gap3_TQtag eq 'T5');
                    $gap3_name1 = $gap3_qname if ($gap3_TQtag eq 'Q5') or ($gap3_TQtag eq 'Q3');
                    ($rate_PE, $read_num) = &find_aligned_read_connect ($gap5_name1, $gap5_TQtag, $gap5_len, $gap3_name1, $gap3_TQtag, $gap3_len, $match_len1);
                }
                if ((($match == 1) and ($PE_align_tag == 1) and (($match_len1 >= $min_match_length) or (($match_len1 >= $min_match_len_local) and ($identity2 >= $min_identity) and (($read_num < 10) or ($rate_PE >= 0.1))) or (($identity2 >= 97) and ($rate_PE >= 0.2)))) or 
                (($match == 1) and ($PE_align_tag == 0) and ($match_len1 >= $min_match_len_local) and ($identity2 >= $min_identity))){
                    if (($gap5_qname ne '') and ($gap3_qname ne '')){
                        if ($match_len1 >= $gap_5_unmatch_len + $gap_3_unmatch_len){
                            my $gap_size = $target_gap_size{$tname};
                            my $overlen = $match_len1 - $gap_5_unmatch_len - $gap_3_unmatch_len;
                            $gap_Cseq{$tname} = $overlen;
                        }
                        elsif (($match_len1 <= $gap_5_unmatch_len) and ($match_len2 <= $gap_3_unmatch_len)){
                            my $gap_seq_1 = substr ($gap5_seq, -$gap_5_unmatch_len, $gap_5_unmatch_len);
                            substr ($gap_seq_1, -$match_len1, $match_len1, '');
                            my $gap_seq_2 = substr ($gap3_seq, 0, $gap_3_unmatch_len);
                            $gap_Aseq{$tname} = $gap_seq_1 . $gap_seq_2;
                            $used_longread_TA{$tname} = "$gap5_qname==$gap3_qname";
                        }
                        elsif ($match_len1 <= $gap_5_unmatch_len){
                            my $gap_seq = substr ($gap5_seq, -$gap_5_unmatch_len, $gap_5_unmatch_len);
                            substr ($gap_seq, -($match_len2 - $gap_3_unmatch_len), $match_len2 - $gap_3_unmatch_len, '');
                            $gap_Aseq{$tname} = $gap_seq;
                            $used_longread_TA{$tname} = "$gap5_qname==$gap3_qname";
                        }
                        elsif ($match_len2 <= $gap_3_unmatch_len){
                            my $gap_seq = substr ($gap3_seq, 0, $gap_3_unmatch_len);
                            substr ($gap_seq, 0, $match_len1 - $gap_5_unmatch_len, '');
                            $gap_Aseq{$tname} = $gap_seq;
                            $used_longread_TA{$tname} = "$gap5_qname==$gap3_qname";
                        }
                    }
                    elsif ($gap5_qname ne ''){
                        if ($match_len1 >= $gap_5_unmatch_len){
                            my $gap_size = $target_gap_size{$tname};
                            my $overlen = $match_len1 - $gap_5_unmatch_len;
                            $gap_Cseq{$tname} = $overlen;
                        }
                        else{
                            my $gap_seq = substr ($gap5_seq, -$gap_5_unmatch_len, $gap_5_unmatch_len);
                            substr ($gap_seq, -$match_len1, $match_len1, '');
                            $gap_Aseq{$tname} = $gap_seq;
                            $used_longread_TA{$tname} = $gap5_qname;
                        }
                    }
                    elsif ($gap3_qname ne ''){
                        if ($match_len2 >= $gap_3_unmatch_len){
                            my $gap_size = $target_gap_size{$tname};
                            my $overlen = $match_len2 - $gap_3_unmatch_len;
                            $gap_Cseq{$tname} = $overlen;
                        }
                        else{
                            my $gap_seq = substr ($gap3_seq, 0, $gap_3_unmatch_len);
                            substr ($gap_seq, 0, $match_len2, '');
                            $gap_Aseq{$tname} = $gap_seq;
                            $used_longread_TA{$down_tname} = $gap3_qname;
                        }
                    }
                    delete $gap_5{$tname} if ($gap5_qname ne '');
                    delete $gap_3{$tname} if ($gap3_qname ne '');
                }
                elsif ($match == 0){
                    my $gap_size = $target_gap_size{$tname};
                    if ($gap_5_unmatch_len + $gap_3_unmatch_len > $gap_size * 1.5){
                        delete $gap_5{$tname} if ($gap5_qname ne '');
                        delete $gap_3{$tname} if ($gap3_qname ne '');
                    }
                }
            }
            else{
                if (($connect_subcon == 1) and ($ite_num == 1)){
                    $gap5_seq = $target_subcon_seq{$tname};
                    $gap5_TQtag = 'T3';
                    $gap3_seq = $target_subcon_seq{$down_tname};
                    $gap3_TQtag = 'T5';
                    next if (length $gap5_seq < 20) or (length $gap3_seq < 20);
                    my ($match, $match_len1, $match_len2, $identity2) = &local_align3 (\$gap5_seq, \$gap3_seq);
                    my $rate_PE = 0;
                    my $read_num = 0;
                    if (($match == 1) and ($PE_align_tag == 1)){
                        my $gap3_len = length $gap3_seq;
                        my $gap5_len = length $gap5_seq;
                        ($rate_PE, $read_num) = &find_aligned_read_connect ($tname, $gap5_TQtag, $gap5_len, $down_tname, $gap3_TQtag, $gap3_len, $match_len1);
                    }
                    my $gap_size = $target_gap_size{$tname};
                    if (($match_len1 >= $min_match_length) or (($PE_align_tag == 1) and ($match_len1 >= $min_match_len_local) and ($identity2 >= $min_identity) and (($read_num < 10) or ($rate_PE >= 0.1))) or (($identity2 >= 97) and ($rate_PE >= 0.2))){
                        $gap_Cseq{$tname} = $match_len1;
                    }
                }
            }
        }
    }
}

my $concat_seq = '';
my $pre_base_name = '';
my $pre_num = 1;
my $out_gapclosed = "$ite_dir$out_prefix.gapclosed.fa";
open (OUT, "> $out_gapclosed");
foreach my $tname (sort keys %target_gap_size){         # substitute target gap regions with query and extend target termini
    next if (exists $skip_subcon{$tname});
    my $base_name = $1 if ($tname =~/(\S+)\/(\d+)$/);
    my $num = $2;
    $num =~ s/^0*//;
    $pre_base_name = $base_name if ($pre_base_name eq '');
    if ($pre_base_name ne $base_name){
        my $last_tname = $pre_base_name . '/' . sprintf ("%04d", $pre_num + 1) if (!exists $target_noscaf{$pre_base_name});
        $last_tname = $pre_base_name . '/' . sprintf ("%04d", $pre_num) if (exists $target_noscaf{$pre_base_name});
        $concat_seq .= $target_subcon_seq{$last_tname};
        if (($extend == 1) and ($ite_num == 1)){
            if (exists $term5{$pre_base_name}){
                my ($qname, $qstart, $qend, $extend_len, $strand) = split (/==/, $term5{$pre_base_name});
                my $ext_seq = substr ($query_seq{$qname}, $qstart - 1, $extend_len);
                if ($strand eq 'Minus'){
                    $ext_seq = reverse $ext_seq;
                    $ext_seq =~ tr/ACGT/TGCA/;
                }
                $concat_seq = $ext_seq . $concat_seq;
            }
            if (exists $term3{$pre_base_name}){
                my ($qname, $qstart, $qend, $extend_len, $strand) = split (/==/, $term3{$pre_base_name});
                my $ext_seq = substr ($query_seq{$qname}, $qstart - 1, $extend_len);
                if ($strand eq 'Minus'){
                    $ext_seq = reverse $ext_seq;
                    $ext_seq =~ tr/ACGT/TGCA/;
                }
                $concat_seq .= $ext_seq;
            }
        }
        print OUT ">$pre_base_name\n";
        print OUT $concat_seq, "\n";
        $concat_seq = '';
    }
    if (!exists $target_noscaf{$base_name}){
        my $gap_seq = '';
        $concat_seq .= $target_subcon_seq{$tname};
        if (exists $gap_A{$tname}){
            my ($qname, $qstart, $qend, $extend_len, $strand) = split (/==/, $gap_A{$tname});
            if ($extend_len > 0){
                $gap_seq = substr ($query_seq{$qname}, $qstart - 1, $extend_len);
                if ($strand eq 'Minus'){
                    $gap_seq = reverse $gap_seq;
                    $gap_seq =~ tr/ACGT/TGCA/;
                }
                $used_longread_TA{$tname} = $qname;
            }
            else{
                my $subcon_overlap_len = 0 - $extend_len;
                my $del = substr ($concat_seq, -$subcon_overlap_len, $subcon_overlap_len, '');
                $gap_seq = 'oo';
            }
        }
        elsif (exists $gap_Aseq{$tname}){
            $gap_seq = $gap_Aseq{$tname};
        }
        elsif (exists $gap_Cseq{$tname}){
            $gap_seq = 'nn';
        }
        else{
            my $gap_seq_5 = '';
            my $gap_seq_3 = '';
            if (exists $gap_5{$tname}){
                my ($qname, $qstart, $qend, $extend_len, $strand) = split (/==/, $gap_5{$tname});
                $gap_seq_5 = substr ($query_seq{$qname}, $qstart - 1, $extend_len);
                if ($strand eq 'Minus'){
                    $gap_seq_5 = reverse $gap_seq_5;
                    $gap_seq_5 =~ tr/ACGT/TGCA/;
                }
                $used_longread_T3{$tname} = $qname;
            }
            if (exists $gap_3{$tname}){
                my ($qname, $qstart, $qend, $extend_len, $strand) = split (/==/, $gap_3{$tname});
                $gap_seq_3 = substr ($query_seq{$qname}, $qstart - 1, $extend_len);
                if ($strand eq 'Minus'){
                    $gap_seq_3 = reverse $gap_seq_3;
                    $gap_seq_3 =~ tr/ACGT/TGCA/;
                }
                my $next_tname = $base_name . '/' . sprintf ("%04d", $num + 1);
                $used_longread_T5{$next_tname} = $qname;
            }
            if ((length $gap_seq_5 > $target_gap_size{$tname} * 1.5) and (length $gap_seq_3 <= $target_gap_size{$tname} * 1.5)){
                if (length $gap_seq_3 <= $target_gap_size{$tname} - $min_gap_size){
                    my $gap_size = $target_gap_size{$tname} - length $gap_seq_3;
                    $gap_seq = ('N' x $gap_size) . $gap_seq_3;
                }
                else{
                    $gap_seq = ('N' x $min_gap_size) . $gap_seq_3 if ($min_gap_size > 5);
                    $gap_seq = 'NNNNN' . $gap_seq_3 if ($min_gap_size <= 5);
                }
                delete $used_longread_T3{$tname};
            }
            elsif ((length $gap_seq_5 <= $target_gap_size{$tname} * 1.5) and (length $gap_seq_3 > $target_gap_size{$tname} * 1.5)){
                if (length $gap_seq_5 <= $target_gap_size{$tname} - $min_gap_size){
                    my $gap_size = $target_gap_size{$tname} - length $gap_seq_5;
                    $gap_seq = $gap_seq_5 . ('N' x $gap_size);
                }
                else{
                    $gap_seq = $gap_seq_5 . ('N' x $min_gap_size) if ($min_gap_size > 5);
                    $gap_seq = $gap_seq_5 . 'NNNNN' if ($min_gap_size <= 5);
                }
                my $next_tname = $base_name . '/' . sprintf ("%04d", $num + 1);
                delete $used_longread_T5{$next_tname};
            }
            elsif (length ($gap_seq_5) + length ($gap_seq_3) <= $target_gap_size{$tname} * 1.5){
                if (length ($gap_seq_5) + length ($gap_seq_3) <= $target_gap_size{$tname} - $min_gap_size){
                    my $gap_size = $target_gap_size{$tname} - length ($gap_seq_5) - length ($gap_seq_3);
                    $gap_seq = $gap_seq_5 . ('N' x $gap_size) . $gap_seq_3;
                }
                else{
                    $gap_seq = $gap_seq_5 . ('N' x $min_gap_size) . $gap_seq_3 if ($min_gap_size > 5);
                    $gap_seq = $gap_seq_5 . 'NNNNN' . $gap_seq_3 if ($min_gap_size <= 5);
                }
            }
        }
        if ($gap_seq eq ''){
            $gap_seq = 'N' x $target_gap_size{$tname};
        }
        if ($gap_seq eq 'nn'){
            my $overlap_len = $gap_Cseq{$tname};
            my $del = substr ($concat_seq, -$overlap_len, $overlap_len, '');
        }
        elsif ($gap_seq ne 'oo'){
            $concat_seq .= $gap_seq;
        }
    }
    $pre_base_name = $base_name;
    $pre_num = $num;
}
my $last_tname = $pre_base_name . '/' . sprintf ("%04d", $pre_num + 1) if (!exists $target_noscaf{$pre_base_name});
$last_tname = $pre_base_name . '/' . sprintf ("%04d", $pre_num) if (exists $target_noscaf{$pre_base_name});
$concat_seq .= $target_subcon_seq{$last_tname};

if (($extend == 1) and ($ite_num == 1)){
    if (exists $term5{$pre_base_name}){
        my ($qname, $qstart, $qend, $extend_len, $strand) = split (/==/, $term5{$pre_base_name});
        my $ext_seq = substr ($query_seq{$qname}, $qstart - 1, $extend_len);
        if ($strand eq 'Minus'){
            $ext_seq = reverse $ext_seq;
            $ext_seq =~ tr/ACGT/TGCA/;
        }
        $concat_seq = $ext_seq . $concat_seq;
    }
    if (exists $term3{$pre_base_name}){
        my ($qname, $qstart, $qend, $extend_len, $strand) = split (/==/, $term3{$pre_base_name});
        my $ext_seq = substr ($query_seq{$qname}, $qstart - 1, $extend_len);
        if ($strand eq 'Minus'){
            $ext_seq = reverse $ext_seq;
            $ext_seq =~ tr/ACGT/TGCA/;
        }
        $concat_seq .= $ext_seq;
    }
}

print OUT ">$pre_base_name\n";
print OUT $concat_seq, "\n";
close (OUT);

my $used_read_file = "$ite_dir$out_prefix-4m.longreads-used-for-gapfill.txt";
open (OUT, ">$used_read_file");
print OUT '-' x 40, "\n";
print OUT "target_name\tregion\tquery_name\n";
print OUT '-' x 40, "\n";
foreach my $tname (sort keys %target_subcon_seq){
    if (exists $used_longread_T5{$tname}){
        print OUT "$tname\t5-term\t$used_longread_T5{$tname}\n";
    }
    if (exists $used_longread_TA{$tname}){
        if ($used_longread_TA{$tname} =~ /==/){
            my @qlist = split (/==/, $used_longread_TA{$tname});
            print OUT "$tname\tgap-all\t@qlist\n";
        }
        else{
            print OUT "$tname\tgap-all\t$used_longread_TA{$tname}\n";
        }
    }
    if (exists $used_longread_T3{$tname}){
        print OUT "$tname\t3-term\t$used_longread_T3{$tname}\n";
    }
}
close (OUT);

my @file_list = <$temp_dir/$out_prefix*>;
system ("rm -f $temp_dir/$out_prefix*") if (@file_list > 0) and ($ite_num == $iteration);

$localtime = localtime;
print STDERR "\nJob at iteration-$ite_num finished: $localtime\n";

################################################################################################################

sub assign_Tsubcont{
    my ($name, $seq) = @_;
    my $sum_size = 0;
    my $count_subcon = 0;
    my $subseq = '';
    my $carry_gap_len = 0;
    $seq = uc $seq;
    $seq =~ s/^N+// if ($seq =~ /^N+/);
    $seq =~ s/N+$// if ($seq =~ /N+$/);
    while ($seq =~ /([A-Z]+?)(N{$min_gap_size,})/g){
        my $subcontig_size = length $1;
        my $gap_size = length $2;
        $subseq = $1 . $2;
        $sum_size += $subcontig_size + $gap_size;
        if (($subcontig_size < $min_subcon_size) and ($count_subcon > 0)){
            $carry_gap_len += $subcontig_size + $gap_size;
            next;
        }
        $count_subcon ++;
        my $count_subcon_04d = sprintf ("%04d", $count_subcon);
        my $subcon_name = $name . '/' . $count_subcon_04d;
        $target_subcon_seq{$subcon_name} = $1;
        $total_target_subcon_len += length $1;
        $target_gap_size{$subcon_name} = $gap_size + $carry_gap_len;
        if (exists $target_Is{$name}){
            $target_Is{$name} .= '=' . ($sum_size - $subcontig_size - $gap_size + 1);
        }
        else{
            $target_Is{$name} = 1;
        }
        if (exists $target_Ie{$name}){
            $target_Ie{$name} .= '=' . ($sum_size - $gap_size);
        }
        else{
            $target_Ie{$name} = $subcontig_size;
        }
        $carry_gap_len = 0;
    }
    if ($seq =~ /$subseq([A-Z]+?)$/){
        $count_subcon ++;
        my $count_subcon_04d = sprintf ("%04d", $count_subcon);
        my $subcon_name = $name . '/' . $count_subcon_04d;
        $target_subcon_seq{$subcon_name} = $1;
        $total_target_subcon_len += length $1;
        $target_Is{$name} .= '=' . ($sum_size + 1);
        $target_Ie{$name} .= '=' . (length $seq);
    }
}

sub filter_blast_result{                        # filter blast results
    my %target_pos;
    foreach my $targetname (keys %target_info){
        foreach my $queryname (keys %{$target_info{$targetname}}){
            my $align = ${$target_info{$targetname}}{$queryname};
            my ($qname, $qstart, $qend, $tname, $tstart, $tend, $strand, $overlen, $identity) = split (/==/, $align);
            my $Qlen = length $query_seq{$qname};
            my $Tlen = length $target_subcon_seq{$tname};
            my $tname_base = $1 if ($tname =~/(\S+)\/(\d+)$/);
            my $tname_num = $2;
            $tname_num =~ s/^0*//;
            if (exists $included_qname{$qname}){
                ${$target_pos{$tname}}{$tstart} = $align . '==' . $Qlen . '==' . $Tlen;
                next;
            }
            if (($tstart <= 3) and ($tend + 2 >= $Tlen)){
                ${$T_match_A{$tname}}{$qname} = $align;
                $included_tname{$tname} = $qname;
                push @{${$TQ_cover{$tname_base}}{$qname}}, $tname_num;
                ${$TQ_strand{$tname_base}}{$qname} .= $strand . '=';
                my $qunmatch_len = $qstart - 1 + $Qlen - $qend;
                ${$Q_match{$qname}}{$tname_base} += $qunmatch_len;
            }
            elsif (($qstart <= 3) and ($qend + 2 >= $Qlen)){
                ${$T_match_I{$tname}}{$qname} = $align;
                $included_qname{$qname} = $tname;
                ${$Q_match{$qname}}{$tname_base} += 0;
            }
            else{
                if ($strand eq 'Plus'){     
                    if (($tstart <= 3) and ($qend + 2 >= $Qlen)){
                        my $tunmatch_len = $Tlen - $tend;
                        my $qunmatch_len = $qstart - 1;
                        if (($PE_align_tag == 1) and (($tunmatch_len > $read_length + 10) or ($qunmatch_len > $read_length + 10))){
                            my ($rate_PE, $read_num) = &find_aligned_read ($align, 'T5');
                            if (($rate_PE > 0) or ($read_num > 5)){
                                $evaluate_pairs{$align} = $rate_PE;
                            }
                            else{
                                $evaluate_pairs{$align} = -1;
                            }
                        }
                        else{
                            $evaluate_pairs{$align} = -1;
                        }
                    }
                    elsif (($tend + 2 >= $Tlen) and ($qstart <= 3)){
                        my $tunmatch_len = $tstart - 1;
                        my $qunmatch_len = $Qlen - $qend;
                        if (($PE_align_tag == 1) and (($tunmatch_len > $read_length + 10) or ($qunmatch_len > $read_length + 10))){
                            my ($rate_PE, $read_num) = &find_aligned_read ($align, 'T3');
                            if (($rate_PE > 0) or ($read_num > 5)){
                                $evaluate_pairs{$align} = $rate_PE;
                            }
                            else{
                                $evaluate_pairs{$align} = -1;
                            }
                        }
                        else{
                            $evaluate_pairs{$align} = -1;
                        }
                    }
                }
                elsif ($strand eq 'Minus'){
                    if (($tstart <= 3) and ($qstart <= 3)){
                        my $tunmatch_len = $Tlen - $tend;
                        my $qunmatch_len = $Qlen - $qend;
                        if (($PE_align_tag == 1) and (($tunmatch_len > $read_length + 10) or ($qunmatch_len > $read_length + 10))){
                            my ($rate_PE, $read_num) = &find_aligned_read ($align, 'T5');
                            if (($rate_PE > 0) or ($read_num > 5)){
                                $evaluate_pairs{$align} = $rate_PE;
                            }
                            else{
                                $evaluate_pairs{$align} = -1;
                            }
                        }
                        else{
                            $evaluate_pairs{$align} = -1;
                        }
                    }
                    elsif (($tend + 2 >= $Tlen) and ($qend + 2 >= $Qlen)){
                        my $tunmatch_len = $tstart - 1;
                        my $qunmatch_len = $qstart - 1;
                        if (($PE_align_tag == 1) and (($tunmatch_len > $read_length + 10) or ($qunmatch_len > $read_length + 10))){
                            my ($rate_PE, $read_num) = &find_aligned_read ($align, 'T3');
                            if (($rate_PE > 0) or ($read_num > 5)){
                                $evaluate_pairs{$align} = $rate_PE;
                            }
                            else{
                                $evaluate_pairs{$align} = -1;
                            }
                        }
                        else{
                            $evaluate_pairs{$align} = -1;
                        }
                    }
                }
            }
            ${$target_pos{$tname}}{$tstart} = $align . '==' . $Qlen . '==' . $Tlen;
        }
    }
    my @output;
    foreach my $tname (sort keys %target_pos){  # sort according query coordinate, remove alignment with identical start or end points
        foreach my $tpos (sort {$a <=> $b} keys %{$target_pos{$tname}}){
            my ($qname, $qstart, $qend, $tname, $tstart, $tend, $strand, $overlen, $identity, $qlen, $tlen) = split (/==/, ${$target_pos{$tname}}{$tpos});
            my $info = $qname . "\t" . $qstart . "\t" . $qend . "\t" . $qlen . "\t" . $tname . "\t" . $tstart . "\t" . $tend . "\t" . $tlen . "\t" . $strand . "\t" . $identity . "%";
            if ((($qstart <= 3) and ($qend >= $qlen - 2)) or (($tstart <= 3) and ($tend >= $tlen - 2))){
                push @output, $info;
            }
            elsif ((($qstart <= 3) and ($tend >=  $tlen - 2) and ($strand eq 'Plus')) or (($tstart <= 3) and ($qend >= $qlen - 2) and ($strand eq 'Plus')) or (($qstart <= 3) and ($tstart <= 3) and ($strand eq 'Minus')) or (($qend >= $qlen - 2) and ($tend >=  $tlen - 2) and ($strand eq 'Minus'))){
                push @output, $info;
            }
            next if ($tstart == 0);
        }
    }
    open (OUT1, ">> $out_blast_filt");
    print OUT1 '-' x 110, "\n";
    print OUT1 "[Query]\t\t\t\t\t[Target]\n";
    print OUT1 "query_name\tpos1\tpos2\tlength\ttarget_name\tpos1\tpos2\tlength\tstrand\tident\n";
    print OUT1 '-' x 110, "\n";
    foreach (@output){
        print OUT1 "$_\n";
    }
    close (OUT1);
    undef @output;
    undef %target_pos;
}

sub filter_blast_result_Q{                        # filter blast results
    my %target_pos;
    foreach my $targetname (keys %target_info){
        foreach my $queryname (keys %{$target_info{$targetname}}){
            my $align = ${$target_info{$targetname}}{$queryname};
            my ($qname, $qstart, $qend, $tname, $tstart, $tend, $strand, $overlen, $identity) = split (/==/, $align);
            my $Qlen = length $query_seq{$qname};
            my $Tlen = length $query_seq{$tname};
            if (exists $included_qname_Q{$qname}){
                ${$target_pos{$tname}}{$tstart} = $align . '==' . $Qlen . '==' . $Tlen;
                next;
            }
            if (($tstart <= 3) and ($tend + 2 >= $Tlen)){
                ${$included_tname_Q{$tname}}{$qname} = $strand;
            }
            elsif (($qstart <= 3) and ($qend + 2 >= $Qlen)){
                ${$included_qname_Q{$qname}}{$tname} = $strand;
            }
            else{
                if ($strand eq 'Plus'){     
                    if (($tstart <= 3) and ($qend + 2 >= $Qlen)){
                        ${$Q_match_5{$tname}}{$qname} = $strand;
                    }
                    elsif (($tend + 2 >= $Tlen) and ($qstart <= 3)){
                        ${$Q_match_3{$tname}}{$qname} = $strand;
                    }
                }
                elsif ($strand eq 'Minus'){
                    if (($tstart <= 3) and ($qstart <= 3)){
                        ${$Q_match_5{$tname}}{$qname} = $strand;
                    }
                    elsif (($tend + 2 >= $Tlen) and ($qend + 2 >= $Qlen)){
                        ${$Q_match_3{$tname}}{$qname} = $strand;
                    }
                }
            }
            ${$target_pos{$tname}}{$tstart} = $align . '==' . $Qlen . '==' . $Tlen;
        }
    }
    my @output;
    foreach my $tname (sort keys %target_pos){  # sort according query coordinate, remove alignment with identical start or end points
        foreach my $tpos (sort {$a <=> $b} keys %{$target_pos{$tname}}){
            my ($qname, $qstart, $qend, $tname, $tstart, $tend, $strand, $overlen, $identity, $qlen, $tlen) = split (/==/, ${$target_pos{$tname}}{$tpos});
            my $info = $qname . "\t" . $qstart . "\t" . $qend . "\t" . $qlen . "\t" . $tname . "\t" . $tstart . "\t" . $tend . "\t" . $tlen . "\t" . $strand . "\t" . $identity . "%";
            if ((($qstart <= 3) and ($qend >= $qlen - 2)) or (($tstart <= 3) and ($tend >= $tlen - 2))){
                push @output, $info;
            }
            elsif ((($qstart <= 3) and ($tend >=  $tlen - 2) and ($strand eq 'Plus')) or (($tstart <= 3) and ($qend >= $qlen - 2) and ($strand eq 'Plus')) or (($qstart <= 3) and ($tstart <= 3) and ($strand eq 'Minus')) or (($qend >= $qlen - 2) and ($tend >=  $tlen - 2) and ($strand eq 'Minus'))){
                push @output, $info;
            }
            next if ($tstart == 0);
        }
    }
    open (OUT1, ">> $ite_dir$out_prefix-query-query.blast-filt.txt");
    print OUT1 '-' x 110, "\n";
    print OUT1 "[Query]\t\t\t\t\t[Target]\n";
    print OUT1 "query_name\tpos1\tpos2\tlength\ttarget_name\tpos1\tpos2\tlength\tstrand\tident\n";
    print OUT1 '-' x 110, "\n";
    foreach (@output){
        print OUT1 "$_\n";
    }
    close (OUT1);
    undef @output;
    undef %target_pos;
}

sub find_aligned_read{
    my ($align_info, $flag) = @_;
    my ($qname, $qstart, $qend, $tname, $tstart, $tend, $strand, $overlen, $ident) = split (/==/, $align_info);
    my $Tlen = length $target_subcon_seq{$tname};
    my $Qlen = length $query_seq{$qname};
    my ($num_pe_suport1, $num_pe_suport2, $num_read1, $num_read2) = (0, 0, 0, 0);
    my $scan_range = $max_insert_size - $read_length;
    my $qscan_start1;
    my $qscan_end1;
    my $qscan_start2;
    my $qscan_end2;
    my $tscan_start1;
    my $tscan_end1;
    my $tscan_start2;
    my $tscan_end2;
    my @target_read_F1 = ();
    my @target_read_R1 = ();
    my @query_read_F1 = ();
    my @query_read_R1 = ();
    my @target_read_F2 = ();
    my @target_read_R2 = ();
    my @query_read_F2 = ();
    my @query_read_R2 = ();
    if ($flag eq 'T5'){
        $tscan_start1 = 1;
        $tscan_end1 = $scan_range - $read_length + 1;
        $tscan_end1 = $Tlen if ($tscan_end1 > $Tlen);
        if ($strand eq 'Plus'){
            $qscan_start1 = $qstart - $scan_range + 1;
            $qscan_start1 = 1 if ($qscan_start1 < 1);
            $qscan_end1 = $qstart - $read_length + 1;
            $qscan_end1 = 1 if ($qscan_end1 < 1);
        }
        elsif ($strand eq 'Minus'){
            $qscan_start1 = $qend + 1;
            $qscan_end1 = $qend + $scan_range - $read_length;
            $qscan_end1 = $Qlen if ($qscan_end1 > $Qlen);
        }
        if (($qscan_end1 - $qscan_start1 >= 10) and ($tscan_end1 - $tscan_start1 >= 10)){
            my @exist_pair;
            my %found;
            for (my $i = $tscan_start1; $i <= $tscan_end1; $i++){
                if (exists ${$target_read_pos_R{$tname}}{$i}){
                    push @target_read_R1, split (/=/, ${$target_read_pos_R{$tname}}{$i});
                }
            }
            if ($strand eq 'Plus'){
                for (my $i = $qscan_start1; $i <= $qscan_end1; $i++){
                    if (exists ${$query_read_pos_F{$qname}}{$i}){
                        push @query_read_F1, split (/=/, ${$query_read_pos_F{$qname}}{$i});
                    }
                }
                map {$found{$_} = 1} @target_read_R1;
                @exist_pair = grep {$found{$_}} @query_read_F1;
                $num_pe_suport1 = scalar @exist_pair;
            }
            else{
                for (my $i = $qscan_start1; $i <= $qscan_end1; $i++){
                    if (exists ${$query_read_pos_R{$qname}}{$i}){
                        push @query_read_R1, split (/=/, ${$query_read_pos_R{$qname}}{$i});
                    }
                }
                map {$found{$_} = 1} @target_read_R1;
                @exist_pair = grep {$found{$_}} @query_read_R1;
                $num_pe_suport1 = scalar @exist_pair;
            }
        }
        $tscan_start2 = $tend + 1;
        $tscan_end2 = $tend + $scan_range - $read_length;
        $tscan_end2 = $Tlen if ($tscan_end2 > $Tlen);
        if ($strand eq 'Plus'){
            $qscan_start2 = $Qlen - $scan_range;
            $qscan_start2 = 1 if ($qscan_start2 < 1);
            $qscan_end2 = $Qlen - $read_length + 1;
            $qscan_end2 = 1 if ($qscan_end2 < 1);
        }
        elsif ($strand eq 'Minus'){
            $qscan_start2 = 1;
            $qscan_end2 = $scan_range - $read_length + 1;
            $qscan_end2 = $Qlen if ($qscan_end2 > $Qlen);
        }
        if (($qscan_end2 - $qscan_start2 >= 10) and ($tscan_end2 - $tscan_start2 >= 10)){
            my @exist_pair;
            my %found;
            for (my $i = $tscan_start2; $i <= $tscan_end2; $i++){
                if (exists ${$target_read_pos_R{$tname}}{$i}){
                    push @target_read_R2, split (/=/, ${$target_read_pos_R{$tname}}{$i});
                }
            }
            if ($strand eq 'Plus'){
                for (my $i = $qscan_start2; $i <= $qscan_end2; $i++){
                    if (exists ${$query_read_pos_F{$qname}}{$i}){
                        push @query_read_F2, split (/=/, ${$query_read_pos_F{$qname}}{$i});
                    }
                }
                map {$found{$_} = 1} @target_read_R2;
                @exist_pair = grep {$found{$_}} @query_read_F2;
                $num_pe_suport2 = scalar @exist_pair;
            }
            else{
                for (my $i = $qscan_start2; $i <= $qscan_end2; $i++){
                    if (exists ${$query_read_pos_R{$qname}}{$i}){
                        push @query_read_R2, split (/=/, ${$query_read_pos_R{$qname}}{$i});
                    }
                }
                map {$found{$_} = 1} @target_read_R2;
                @exist_pair = grep {$found{$_}} @query_read_R2;
                $num_pe_suport2 = scalar @exist_pair;
            }
        }
    }
    elsif ($flag eq 'T3'){
        $tscan_start1 = $Tlen - $scan_range;
        $tscan_start1 = 1 if ($tscan_start1 < 1);
        $tscan_end1 = $Tlen - $read_length + 1;
        $tscan_end1 = 1 if ($tscan_end1 < 1);
        if ($strand eq 'Plus'){
            $qscan_start1 = $qend + 1;
            $qscan_end1 = $qend + $scan_range - $read_length;
            $qscan_end1 = $Qlen if ($qscan_end1 > $Qlen);
        }
        elsif ($strand eq 'Minus'){
            $qscan_start1 = $qstart - $scan_range + 1;
            $qscan_start1 = 1 if ($qscan_start1 < 1);
            $qscan_end1 = $qstart - $read_length + 1;
            $qscan_end1 = 1 if ($qscan_end1 < 1);
        }
        if (($qscan_end1 - $qscan_start1 >= 10) and ($tscan_end1 - $tscan_start1 >= 10)){
            my @exist_pair;
            my %found;
            for (my $i = $tscan_start1; $i <= $tscan_end1; $i++){
                if (exists ${$target_read_pos_F{$tname}}{$i}){
                    push @target_read_F1, split (/=/, ${$target_read_pos_F{$tname}}{$i});
                }
            }
            if ($strand eq 'Plus'){
                for (my $i = $qscan_start1; $i <= $qscan_end1; $i++){
                    if (exists ${$query_read_pos_R{$qname}}{$i}){
                        push @query_read_R1, split (/=/, ${$query_read_pos_R{$qname}}{$i});
                    }
                }
                map {$found{$_} = 1} @target_read_F1;
                @exist_pair = grep {$found{$_}} @query_read_R1;
                $num_pe_suport1 = scalar @exist_pair;
            }
            else{
                for (my $i = $qscan_start1; $i <= $qscan_end1; $i++){
                    if (exists ${$query_read_pos_F{$qname}}{$i}){
                        push @query_read_F1, split (/=/, ${$query_read_pos_F{$qname}}{$i});
                    }
                }
                map {$found{$_} = 1} @target_read_F1;
                @exist_pair = grep {$found{$_}} @query_read_F1;
                $num_pe_suport1 = scalar @exist_pair;
            }
        }
        $tscan_start2 = $tstart - $scan_range + 1;
        $tscan_start2 = 1 if ($tscan_start2 < 1);
        $tscan_end2 = $tstart - $read_length + 1;
        $tscan_end2 = 1 if ($tscan_end2 < 1);
        if ($strand eq 'Plus'){
            $qscan_start2 = 1;
            $qscan_end2 = $scan_range - $read_length + 1;
            $qscan_end2 = $Qlen if ($qscan_end2 > $Qlen);
        }
        elsif ($strand eq 'Minus'){
            $qscan_start2 = $Qlen - $scan_range;
            $qscan_start2 = 1 if ($qscan_start2 < 1);
            $qscan_end2 = $Qlen - $read_length + 1;
            $qscan_end2 = 1 if ($qscan_end2 < 1);
        }
        if (($qscan_end2 - $qscan_start2 >= 10) and ($tscan_end2 - $tscan_start2 >= 10)){
            my @exist_pair;
            my %found;
            for (my $i = $tscan_start2; $i <= $tscan_end2; $i++){
                if (exists ${$target_read_pos_F{$tname}}{$i}){
                    push @target_read_F2, split (/=/, ${$target_read_pos_F{$tname}}{$i});
                }
            }
            if ($strand eq 'Plus'){
                for (my $i = $qscan_start2; $i <= $qscan_end2; $i++){
                    if (exists ${$query_read_pos_R{$qname}}{$i}){
                        push @query_read_R2, split (/=/, ${$query_read_pos_R{$qname}}{$i});
                    }
                }
                map {$found{$_} = 1} @target_read_F2;
                @exist_pair = grep {$found{$_}} @query_read_R2;
                $num_pe_suport2 = scalar @exist_pair;
            }
            else{
                for (my $i = $qscan_start2; $i <= $qscan_end2; $i++){
                    if (exists ${$query_read_pos_F{$qname}}{$i}){
                        push @query_read_F2, split (/=/, ${$query_read_pos_F{$qname}}{$i});
                    }
                }
                map {$found{$_} = 1} @target_read_F2;
                @exist_pair = grep {$found{$_}} @query_read_F2;
                $num_pe_suport2 = scalar @exist_pair;
            }
        }
    }
    my $tscan_range1 = $tscan_end1 - $tscan_start1 + 1;
    my $tscan_range2 = $tscan_end2 - $tscan_start2 + 1;
    my $qscan_range1 = $qscan_end1 - $qscan_start1 + 1;
    my $qscan_range2 = $qscan_end2 - $qscan_start2 + 1;
    if ($tscan_range1 <= $qscan_range1){
        $num_read1 = @target_read_F1 + @target_read_R1;
    }
    else{
        $num_read1 = @query_read_F1 + @query_read_R1;
    }
    if ($tscan_range2 <= $qscan_range2){
        $num_read2 = @target_read_F2 + @target_read_R2;
    }
    else{
        $num_read2 = @query_read_F2 + @query_read_R2;
    }
    my $total_aligned_read = $num_read1 + $num_read2;
    my $rate_PE = 0;
    $rate_PE = int (($num_pe_suport1 + $num_pe_suport2) / $total_aligned_read * 1000) / 1000 if ($total_aligned_read > 0);
    return ($rate_PE, $total_aligned_read);
}

sub find_aligned_read_connect{
    my ($up_cont, $up_tag, $up_len, $cont, $tag, $len, $overlen) = @_;
    my $up_tcont = '';
    my $tcont = '';
    my ($num_pe_suport1, $num_pe_suport2, $num_read1, $num_read2) = (0, 0, 0, 0);
    my $scan_range = $max_insert_size - $read_length;
    my $up_scan_start1;
    my $up_scan_end1;
    my $up_scan_start2;
    my $up_scan_end2;
    my $scan_start1;
    my $scan_end1;
    my $scan_start2;
    my $scan_end2;
    my @up_read_F1 = ();
    my @up_read_R1 = ();
    my @read_F1 = ();
    my @read_R1 = ();
    my @up_read_F2 = ();
    my @up_read_R2 = ();
    my @read_F2 = ();
    my @read_R2 = ();
    if (($tag eq 'T5') or ($tag eq 'Q5')){
        $scan_start1 = 1;
        $scan_end1 = $scan_range - $read_length + 1;
        $scan_end1 = $len if ($scan_end1 > $len);
    }
    elsif ($tag eq 'Q3'){
        $scan_start1 = $len - $scan_range;
        $scan_start1 = 1 if ($len <= $scan_range);
        $scan_end1 = $len - $read_length + 1;
    }
    if (($up_tag eq 'T3') or ($up_tag eq 'Q3')){
        $up_scan_start1 = $up_len - $overlen - $scan_range + 1;
        $up_scan_start1 = 1 if ($up_scan_start1 < 1);
        $up_scan_end1 = $up_len - $overlen - $read_length + 1;
        $up_scan_end1 = 1 if ($up_scan_end1 < 1);
    }
    elsif ($up_tag eq 'Q5'){
        $up_scan_start1 = $overlen + 1;
        $up_scan_end1 = $overlen + $scan_range;
        $up_scan_end1 = $up_len if ($up_scan_end1 > $up_len);
    }
    
    if (($scan_end1 - $scan_start1 >= 10) and ($up_scan_end1 - $up_scan_start1 >= 10)){
        my @exist_pair;
        my %found;
        for (my $i = $up_scan_start1; $i <= $up_scan_end1; $i++){
            if ($up_tag eq 'T3'){
                if (exists ${$target_read_pos_F{$up_cont}}{$i}){
                    push @up_read_F1, split (/=/, ${$target_read_pos_F{$up_cont}}{$i});
                }
            }
            elsif ($up_tag eq 'Q3'){
                if (exists ${$query_read_pos_F{$up_cont}}{$i}){
                    push @up_read_F1, split (/=/, ${$query_read_pos_F{$up_cont}}{$i});
                }
            }
            elsif ($up_tag eq 'Q5'){
                if (exists ${$query_read_pos_R{$up_cont}}{$i}){
                    push @up_read_R1, split (/=/, ${$query_read_pos_R{$up_cont}}{$i});
                }
            }
        }
        for (my $i = $scan_start1; $i <= $scan_end1; $i++){
            if ($tag eq 'T5'){
                if (exists ${$target_read_pos_R{$cont}}{$i}){
                    push @read_R1, split (/=/, ${$target_read_pos_R{$cont}}{$i});
                }
            }
            elsif ($tag eq 'Q5'){
                if (exists ${$query_read_pos_R{$cont}}{$i}){
                    push @read_R1, split (/=/, ${$query_read_pos_R{$cont}}{$i});
                }
            }
            elsif ($tag eq 'Q3'){
                if (exists ${$query_read_pos_F{$cont}}{$i}){
                    push @read_F1, split (/=/, ${$query_read_pos_F{$cont}}{$i});
                }
            }
        }
        my @up_read1 = (@up_read_F1, @up_read_R1);
        my @read1 = (@read_F1, @read_R1);
        map {$found{$_} = 1} @up_read1;
        @exist_pair = grep {$found{$_}} @read1;
        $num_pe_suport1 = scalar @exist_pair;
    }
    
    if (($tag eq 'T5') or ($tag eq 'Q5')){
        $scan_start2 = $overlen + 1;
        $scan_end2 = $scan_range + $overlen - $read_length + 2;
        $scan_end2 = $len if ($scan_end1 > $len) and ($tag eq 'Q5');
        $scan_end2 = $len if ($scan_end1 > $len) and ($tag eq 'T5');
    }
    elsif ($tag eq 'Q3'){
        $scan_start2 = $len - $overlen - $scan_range + 1;
        $scan_start2 = 1 if ($scan_start1 < 1);
        $scan_end2 = $len - $overlen - $read_length + 1;
        $scan_end2 = 1 if ($scan_end1 < 1);
    }
    if (($up_tag eq 'T3') or ($up_tag eq 'Q3')){
        $up_scan_start2 = $up_len - $scan_range + 1 if ($up_tag eq 'Q3');
        $up_scan_start2 = $up_len - $scan_range + 1 if ($up_tag eq 'T3');
        $up_scan_start2 = 1 if ($up_scan_start1 < 1);
        $up_scan_end2 = $up_len - $read_length + 1 if ($up_tag eq 'Q3');
        $up_scan_end2 = $up_len - $read_length + 1 if ($up_tag eq 'T3');
        $up_scan_end2 = 1 if ($up_scan_end1 < 1);
    }
    elsif ($up_tag eq 'Q5'){
        $up_scan_start2 = 1;
        $up_scan_end2 = $scan_range;
        $up_scan_end2 = $up_len if ($scan_range > $up_len);
    }
    
    if (($scan_end2 - $scan_start2 >= 10) and ($up_scan_end2 - $up_scan_start2 >= 10)){
        my @exist_pair = ();
        my %found = ();
        for (my $i = $up_scan_start2; $i <= $up_scan_end2; $i++){
            if ($up_tag eq 'T3'){
                if (exists ${$target_read_pos_F{$up_cont}}{$i}){
                    push @up_read_F2, split (/=/, ${$target_read_pos_F{$up_cont}}{$i});
                }
            }
            elsif ($up_tag eq 'Q3'){
                if (exists ${$query_read_pos_F{$up_cont}}{$i}){
                    push @up_read_F2, split (/=/, ${$query_read_pos_F{$up_cont}}{$i});
                }
            }
            elsif ($up_tag eq 'Q5'){
                if (exists ${$query_read_pos_R{$up_cont}}{$i}){
                    push @up_read_R2, split (/=/, ${$query_read_pos_R{$up_cont}}{$i});
                }
            }
        }
        for (my $i = $scan_start2; $i <= $scan_end2; $i++){
            if ($tag eq 'T5'){
                if (exists ${$target_read_pos_R{$cont}}{$i}){
                    push @read_R2, split (/=/, ${$target_read_pos_R{$cont}}{$i});
                }
            }
            elsif ($tag eq 'Q5'){
                if (exists ${$query_read_pos_R{$cont}}{$i}){
                    push @read_R2, split (/=/, ${$query_read_pos_R{$cont}}{$i});
                }
            }
            elsif ($tag eq 'Q3'){
                if (exists ${$query_read_pos_F{$cont}}{$i}){
                    push @read_F2, split (/=/, ${$query_read_pos_F{$cont}}{$i});
                }
            }
        }
        my @up_read2 = (@up_read_F2, @up_read_R2);
        my @read2 = (@read_F2, @read_R2);
        map {$found{$_} = 1} @up_read2;
        @exist_pair = grep {$found{$_}} @read2;
        $num_pe_suport2 = scalar @exist_pair;
    }
    
    my $up_scan_range1 = $up_scan_end1 - $up_scan_start1 + 1;
    my $up_scan_range2 = $up_scan_end2 - $up_scan_start2 + 1;
    my $scan_range1 = $scan_end1 - $scan_start1 + 1;
    my $scan_range2 = $scan_end2 - $scan_start2 + 1;
    if ($up_scan_range1 <= $scan_range1){
        $num_read1 = @up_read_F1 + @up_read_R1;
    }
    else{
        $num_read1 = @read_F1 + @read_R1;
    }
    if ($up_scan_range2 <= $scan_range2){
        $num_read2 = @up_read_F2 + @up_read_R2;
    }
    else{
        $num_read2 = @read_F2 + @read_R2;
    }
    my $total_aligned_read = $num_read1 + $num_read2;
    my $rate_PE = 0;
    $rate_PE = int (($num_pe_suport1 + $num_pe_suport2) / $total_aligned_read * 1000) / 1000 if ($total_aligned_read > 0);
    return ($rate_PE, $total_aligned_read);
}

sub validate_align{
    my ($overlen, $identity, $rate_PE) = @_;
    my $score = 0;
    if ($overlen < 70){
        $score = 0.052;
    }
    elsif ($overlen < 100){
        $score = 0.076;
    }
    elsif ($overlen < 150){
        $score = 0.17;
    }
    elsif ($overlen < 200){
        $score = 0.4;
    }
    elsif ($overlen < 300){
        $score = 0.93;
    }
    elsif ($overlen < 500){
        $score = 3.1;
    }
    elsif ($overlen >= 500){
        $score = 13.2;
    }
    if ($identity < 96){
        $score *= 0.017;
    }
    elsif ($identity < 97){
        $score *= 0.022;
    }
    elsif ($identity < 98){
        $score *= 0.027;
    }
    elsif ($identity < 99){
        $score *= 0.054;
    }
    elsif ($identity >= 99){
        $score *= 3.08;
    }
    if ($rate_PE < 0.001){
        $score *= 0.095;
    }
    elsif ($rate_PE < 0.01){
        $score *= 0.15;
    }
    elsif ($rate_PE < 0.02){
        $score *= 0.36;
    }
    elsif ($rate_PE < 0.05){
        $score *= 1.03;
    }
    elsif ($rate_PE < 0.1){
        $score *= 4.1;
    }
    elsif ($rate_PE < 0.2){
        $score *= 15.5;
    }
    elsif ($rate_PE >= 0.2){
        $score *= 43;
    }
    my $log_score = int (log ($score) * 100) / 100;
    if (($overlen >= 0) and ($identity >= 0)){
        $num_score ++;
        $sum_score += $score;
        if ($rate_PE > 0){
            $num_score_2 ++;
            $sum_score_2 += $score;
        }
    }
    return ($log_score);
}

sub connect_subcon{
    my $subfile = shift;
    my $num = shift;
    my $out_prefix_2 = shift;
    my $temp_dir_2 = shift;
    my $queryfile = shift;
    my $gapA_file = shift;
    my $gap5_file = shift;
    my $gap3_file = shift;
    my $pe_align_tag = shift;
    my $max_insert_size_2 = shift if ($pe_align_tag == 1);
    my $read_length_2 = shift if ($pe_align_tag == 1);
    my $T_PEalignF_file = shift if ($pe_align_tag == 1);
    my $T_PEalignR_file = shift if ($pe_align_tag == 1);
    my $Q_PEalignF_file = shift if ($pe_align_tag == 1);
    my $Q_PEalignR_file = shift if ($pe_align_tag == 1); 
    my $connect_file = "$temp_dir/$out_prefix-subcon-connect-$num.txt";
    my $delete_tqname_file = "$temp_dir/$out_prefix-delete_tqname-$num.txt";
    if ($connect_subcon == 0){
        system ("$Bin/connect_subcontigs_GMcloser2.pl -s $subfile -q $queryfile -n $num -p $out_prefix_2 -td $temp_dir_2 -ga $gapA_file -gu $gap5_file -gd $gap3_file -t $pe_align_tag -mi $min_identity -mm $min_match_length -ml $min_match_len_local -in $ite_num -i $max_insert_size_2 -r $read_length_2 -tf $T_PEalignF_file -tr $T_PEalignR_file -qf $Q_PEalignF_file -qr $Q_PEalignR_file") if ($pe_align_tag == 1);
        system ("$Bin/connect_subcontigs_GMcloser2.pl -s $subfile -q $queryfile -n $num -p $out_prefix_2 -td $temp_dir_2 -ga $gapA_file -gu $gap5_file -gd $gap3_file -t $pe_align_tag -mi $min_identity -mm $min_match_length -ml $min_match_len_local -in $ite_num") if ($pe_align_tag == 0);
    }
    else{
        system ("$Bin/connect_subcontigs_GMcloser2.pl -s $subfile -q $queryfile -n $num -p $out_prefix_2 -td $temp_dir_2 -ga $gapA_file -gu $gap5_file -gd $gap3_file -t $pe_align_tag -mi $min_identity -mm $min_match_length -ml $min_match_len_local -in $ite_num -i $max_insert_size_2 -r $read_length_2 -tf $T_PEalignF_file -tr $T_PEalignR_file -qf $Q_PEalignF_file -qr $Q_PEalignR_file -cn") if ($pe_align_tag == 1);
        system ("$Bin/connect_subcontigs_GMcloser2.pl -s $subfile -q $queryfile -n $num -p $out_prefix_2 -td $temp_dir_2 -ga $gapA_file -gu $gap5_file -gd $gap3_file -t $pe_align_tag -mi $min_identity -mm $min_match_length -ml $min_match_len_local -in $ite_num -cn") if ($pe_align_tag == 0);
    }
    threads->yield();
    sleep 1;
    return ($connect_file, $delete_tqname_file, $num);
}

sub local_align3{                    # align between two query sequences
    my ($ref_seq1, $ref_seq2) = @_;        # $qseq1: upstream subcon + unamtch region of the overlapped query, $qseq2: unmatch region of query that overlaps target subcon
    my $seq1 = ${$ref_seq1};
    my $seq2 = ${$ref_seq2};
    open (NEWFILE1, "> temp/$out_prefix-seq1.fa");
    print NEWFILE1 '>seq1', "\n", $seq1;
    close (NEWFILE1);
                
    open (NEWFILE1, "> temp/$out_prefix-seq2.fa");
    print NEWFILE1 '>seq2', "\n", $seq2;
    close (NEWFILE1);
    my $mutation_rate = 3;
    my @result = `yass -O 5 -m $mutation_rate -r 0 temp/$out_prefix-seq1.fa temp/$out_prefix-seq2.fa 2>/dev/null`;
    my $match_len1 = 0;
    my $match_len2 = 0;
    my $matchbase_num = 0;
    my $match = 0;
    foreach my $line (@result){
        chomp $line;
        if (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match == 0)){
            $match_len1 = $2 - $1 + 1;
            $match_len2 = $4 - $3 + 1;
            if (($2 >= length ($seq1) - 1) and ($3 <= 2)){
                $match = 1;
                if ($2 == length ($seq1) - 1){
                    $match_len1 ++;
                    $match_len2 ++;
                }
                if ($3 == 2){
                    $match_len1 ++;
                    $match_len2 ++;
                }
            }
        }
        elsif (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match == 1)){
            last;
        }
        if ($match == 1){
            if ($line =~ /^[\s\.\|:]+$/){
                my $count_match = $line =~ s/\|/\|/g;
                $matchbase_num += $count_match;
            }
        }
    }
    undef @result;
    if (($matchbase_num < 20) or ($match_len1 < 20)){
        return (0, 0, 0, 0);
    }
    elsif ($match == 1){
        my $identity = int ($matchbase_num / $match_len1 * 100);
        return ($match, $match_len1, $match_len2, $identity);
    }
    else{
        return (0, 0, 0, 0);
    }
}
