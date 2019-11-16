#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use Config;
use threads;

# a script for consecutive commands of Nucmer/GAGE and correct_contigs_GAGE_coords.pl for evaluation of contig integrity accuracy

my $contig_file = '';
my $ref_file = '';
my $min_identity = 97;
my $min_contig_len = 200;
my $min_coverage = 99;
my $prefix_out = 'out';
my $error_correct = 0;
my $min_nucmer_match = 30;
my $max_indel_len = 100;
my $thread = 1;
my $help;

GetOptions(
    'query|q=s' => \$contig_file,
    'ref|r=s' => \$ref_file,
    'min_id|mi=i' => \$min_identity,
    'min_len|ml=i' => \$min_contig_len,
    'min_cov|mc=i' => \$min_coverage,
    'prefix|p=s' => \$prefix_out,
    'enable_out|e' => \$error_correct,
    'nuc_len|l=i' => \$min_nucmer_match,
    'max_indel|is=i' => \$max_indel_len,
    'thread|n=i' => \$thread,
    'help' => \$help
) or pod2usage(-verbose => 0);

pod2usage(-verbose => 0) if $help;
die "--query option is not specifies: $!\n" if ($contig_file eq '');
die "--ref option is not specified: $!\n" if ($ref_file eq '');
#die "The string specified with --prefix option must not contain directory names but only a file prefix name: $!\n" if ($prefix_out =~ /\//);

=head1 SYNOPSIS
  
  Usage: gmvalid contig -r [reference.fasta] -q [contig.fasta] -p [prefix name of output] [other options]
  Options:
   --query or -q <STR>      input contig fasta file (e.g., contig1.fa)
   --ref or -r <STR>        input reference file (e.g., ref.fa)
   --min_id or -mi <INT>    minimum alignment identity (%) [default: 97]
   --min_len or -ml <INT>   minimum contig length to be considered [default: 200]
   --min_cov or -mc <INT>   minimum coverage (%) of query (contig) aligned to a reference [default: 99]
   --prefix or -p <STR>     prefix name of output files
   --error_correct or -e    output an error-corrected contig set [default: false]
   --nuc_len or -l <INT>    minimum exact match length for specifying nucmer option -l [default: 30]
   --max_indel or -is <INT> maximum allowable size of indels (or distance between break points of a local misassembly) [default: 100]
   --thread or -n           number of threads to run [default: 1]
   --help or -h             output help message
   
=cut


my $temp_dir = 'temp';
my $prefix_out_base = $prefix_out;
my $out_dir = '';
if ($prefix_out =~ /^(\S+)\/([^\/]+)$/){
    $out_dir = $1;
    $prefix_out_base = $2;
    $temp_dir = "$out_dir/temp";
    system ("mkdir $out_dir") unless (-d $out_dir);
}
system ("mkdir $temp_dir") unless (-d $temp_dir);

$thread = 24 if ($thread > 24);
$thread = 1 if (!$Config{useithreads});

my @target;
my $seq = '';
my $target_length = 0;

open (FILE, $ref_file) or die "$!";
while (my $line = <FILE>){
    if ($line =~ /^>/){
        if ($seq ne ''){
            $target_length += length $seq;
            push @target, $seq;
        }
        $seq = $line;
    }
    else{
        $seq .= $line;
    }
}
$target_length += length $seq;
push @target, $seq;
$seq = '';
close (FILE);

if ($min_nucmer_match == 30){
    $min_nucmer_match = 40 if ($target_length >= 200000000);
    $min_nucmer_match = 50 if ($target_length >= 300000000);
}

my $max_subseq_len = 50000000;
my $sec_num = 1;
my $dev_len = 0;
while (1){
    if ($target_length / $sec_num / $thread <= $max_subseq_len){
        $dev_len = int ($target_length / $sec_num / $thread) + $thread;
        last;
    }
    else{
        $sec_num ++;
    }
}

my $count_sub = 1;
my $sum_subseq = 0;
my @target_subfiles;
my @sub_target;

if (@target > $thread){
    foreach my $chr (@target){
        push @sub_target, $chr;
        $sum_subseq += length $chr;
        if ($sum_subseq >= $dev_len){
            open (OUT, "> $temp_dir/$prefix_out_base-target-$count_sub.fasta");
            foreach (@sub_target){
                print OUT $_;
            }
            close (OUT);
            @sub_target = ();
            push @target_subfiles, "$temp_dir/$prefix_out_base-target-$count_sub.fasta";
            $sum_subseq = 0;
            $count_sub ++;
        }
    }
    if (@sub_target > 0){
        if ($count_sub > 1){
            open (OUT, "> $temp_dir/$prefix_out_base-target-$count_sub.fasta");
            foreach (@sub_target){
                print OUT $_;
            }
            close (OUT);
            push @target_subfiles, "$temp_dir/$prefix_out_base-target-$count_sub.fasta";
        }
        else{
            push @target_subfiles, $ref_file;
        }
    }
}
else{
    foreach my $chr (@target){
        push @sub_target, $chr;
        open (OUT, "> $temp_dir/$prefix_out_base-target-$count_sub.fasta");
        foreach (@sub_target){
            print OUT $_;
        }
        close (OUT);
        @sub_target = ();
        push @target_subfiles, "$temp_dir/$prefix_out_base-target-$count_sub.fasta";
        $count_sub ++;
    }
}
undef @sub_target;
undef @target;

my @nucmer_subfiles;

if ($thread > 1){
    my @align_jobs;
    my $count = 0;
    foreach my $subfile (@target_subfiles){
        $count ++;
        my ($thread_job) = threads->new(\&nucmer_align, $subfile);
        push (@align_jobs, $thread_job);
        if ($count == $thread){
            foreach (@align_jobs){
                my ($file_name) = $_->join;
                print STDERR "Nucmer completed ($file_name)\n";
                push @nucmer_subfiles, $file_name;
            }
            @align_jobs = ();
            $count = 0;
        }
    }
    undef @target_subfiles;
    if (@align_jobs > 0){
        foreach (@align_jobs){
            my ($file_name) = $_->join;
            print STDERR "Nucmer completed ($file_name)\n";
            push @nucmer_subfiles, $file_name;
        }
    }
}
else{
    if (@target_subfiles == 1){
        system ("nucmer --maxmatch -p $prefix_out -l $min_nucmer_match -banded -D 5 $ref_file $contig_file");
    }
    elsif (@target_subfiles > 1){
        foreach my $subfile (@target_subfiles){
            my $file_prefix = $1 if ($subfile =~ /(\S+)\.fasta$/);
            system ("nucmer --maxmatch -p $file_prefix -l $min_nucmer_match -banded -D 5 $subfile $contig_file");
            push @nucmer_subfiles, "$file_prefix.delta";
            print STDERR "Nucmer completed ($file_prefix.delta)\n";
        }
    }
}

if (@nucmer_subfiles > 0){
    open (OUT, "> $prefix_out.delta");
    print OUT "$ref_file $contig_file\n";
    print OUT "NUCMER\n";
    
    foreach my $outfile (@nucmer_subfiles){
        open (FILE, $outfile) or die "$outfile is not found: $!\n";
        my $count = 0;
        while (my $line = <FILE>){
            $count ++;
            chomp $line;
            next if ($count <= 2) or ($line =~ /^$/);
            print OUT "$line\n";
        }
        close (FILE);
    }
    close (OUT);
    system ("rm -f $temp_dir/$prefix_out_base-target*");
}

system ("delta-filter -o 95 -i $min_identity $prefix_out.delta > $prefix_out.fdelta");
system ("show-coords -THrcl $prefix_out.fdelta > $prefix_out.coords");

sub nucmer_align{
    my $subfile = shift;
    my $file_prefix = $1 if ($subfile =~ /(\S+)\.fasta$/);
    system ("nucmer --maxmatch -p $file_prefix -l $min_nucmer_match -banded -D 5 $subfile $contig_file 2> /dev/null");
    threads->yield();
#    sleep 5;
    return ("$file_prefix.delta");
}
