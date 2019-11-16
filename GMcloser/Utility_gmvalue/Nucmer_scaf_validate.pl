#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use Config;
use threads;

# a script for consecutive commands of Nucmer/GAGE and correct_scafs_GAGE_coords.pl for evaluation of scaffold contiguity accuracy

my $scaf_file = '';
my $ref_file = '';
my $min_identity = 97;
my $min_coverage = 99;
my $min_subcon_len = 200;
my $min_align_len = 200;
my $min_gap_len = 1;
my $max_gap_len = 50000;
my $max_indel_size = 100;
my $prefix_out = 'out';
my $error_correct = 0;
my $min_nucmer_match = 30;
my $thread = 1;
my $help;

GetOptions(
    'query|q=s' => \$scaf_file,
    'ref|r=s' => \$ref_file,
    'min_id|mi=i' => \$min_identity,
    'min_cov|mc=i' => \$min_coverage,
    'min_len|ml=i' => \$min_subcon_len,
    'min_align|ma=i' => \$min_align_len,
    'min_gap|g=i' => \$min_gap_len,
    'max_gap|mg=i' => \$max_gap_len,
    'max_indel|is=i' => \$max_indel_size,
    'prefix|p=s' => \$prefix_out,
    'enable_out|e' => \$error_correct,
    'nuc_len|l' => \$min_nucmer_match,
    'thread|n=i' => \$thread,
    'help' => \$help
) or pod2usage(-verbose => 0);

pod2usage(-verbose => 0) if $help;
die "--query option is not specifies: $!\n" if ($scaf_file eq '');
die "--ref option is not specified: $!\n" if ($ref_file eq '');

=head1 SYNOPSIS

  GMvalue ver. 1.3

  Usage: gmvalue scaf -r [reference.fasta] -q [scaffold.fasta] -p [prefix name of output] [other options]
  Options:
   --query or -q <STR>      input scaffold fasta file (e.g., scaf1.fa)
   --ref or -r <STR>        input reference file (e.g., ref.fa)
   --min_id or -mi <INT>    minimum alignment identity (%) [default: 97]
   --min_cov or -mc <INT>   minimum coverage (%) of query (contig) aligned to a reference [default: 99]
   --min_align or -ma <INT> minimum alignment overlap length with the maximum allowable size of indels [default: 200]
   --min_len or -ml <INT>   minimum contig length to be considered [default: 200]
   --prefix or -p <STR>     prefix name of output files
   --error_correct or -e    output an error-corrected contig set [default: false]
   --nuc_len or -l <INT>    minimum exact match length for specifying nucmer option -l [default: 30]
   --min_gap or -mg <INT>   minimum gap size in query scaffolds to split into subcontigs [default: 1]
   --max_gap or -mg <INT>   maximum length of gaps contained in the scaffolds [default: 50000]
   --max_indel or -is <INT> maximum allowable size of indels in subcontigs (or distance between break points of a local misassembly) [default: 100]
   --thread or -n           number of threads to run [default: 1]
   --help or -h             output help message
   
=cut

#my $subcontig_file = "$prefix_out.scaf2cont.fasta";
#system ("split_scaf2contig.pl -g $min_gap_len -q $scaf_file > $subcontig_file");

my $subcontig_file = "$prefix_out.subcontig.fasta";
#system ("java -cp ~/tools/gage-paper-validation SplitFastaByLetter $scaf_file N > $subcontig_file");

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

my $seq_name = '';
my $seq = '';
my $query_length = 0;

open (FILE, $scaf_file) or die "$scaf_file is not found: $!\n";
open (OUT, "> $subcontig_file");
while (my $line = <FILE>){
    $line =~ s/(\r\n|\n|\r)//g;
    if ($line =~ /^>/){
        if ($seq ne ''){
            $seq =~ s/^[NnX]+// if ($seq =~ /^[NnX]/);
            $seq =~ s/[NnX]+$// if ($seq =~ /[NnX]$/);
            if ($seq =~ /[Nn]{$min_gap_len}/){
                my $subseq = '';
                my $sum_size = 0;
                my $count_subcon = 0;
                while ($seq =~ /([A-Z]+?)(N{$min_gap_len,})/g){
                    my $subcontig = $1;
                    $subseq = $1 . $2;
                    my $subcon_name = "$seq_name/$count_subcon";
                    $count_subcon ++;
                    $query_length += length $subcontig;
                    print OUT ">$subcon_name\n";
                    print OUT $subcontig, "\n";
                }
                if ($seq =~ /$subseq([A-Z]+?)$/){
                    my $subcon_name = "$seq_name/$count_subcon";
                    my $subcontig = $1;
                    $query_length += length $subcontig;
                    print OUT ">$subcon_name\n";
                    print OUT $subcontig, "\n";
                }
            }
            else{
                my $subcon_name = "$seq_name/0";
                $query_length += length $seq;
                print OUT ">$subcon_name\n";
                print OUT $seq, "\n";
            }
            $seq = '';
        }
        $seq_name = $1 if ($line =~ /^>(\S+)/);
    }
    else{
        $seq .= uc $line;
    }
}
$seq =~ s/^[NnX]+// if ($seq =~ /^[NnX]/);
$seq =~ s/[NnX]+$// if ($seq =~ /[NnX]$/);
if ($seq =~ /[Nn]{$min_gap_len}/){
    my $subseq = '';
    my $sum_size = 0;
    my $count_subcon = 0;
    while ($seq =~ /([A-Z]+?)(N{$min_gap_len,})/g){
        my $subcontig = $1;
        $subseq = $1 . $2;
        my $subcon_name = "$seq_name/$count_subcon";
        $count_subcon ++;
        $query_length += length $subcontig;
        print OUT ">$subcon_name\n";
        print OUT $subcontig, "\n";
    }
    if ($seq =~ /$subseq([A-Z]+?)$/){
        my $subcon_name = "$seq_name/$count_subcon";
        my $subcontig = $1;
        $query_length += length $subcontig;
        print OUT ">$subcon_name\n";
        print OUT $subcontig, "\n";
    }
}
else{
    my $subcon_name = "$seq_name/0";
    $query_length += length $seq;
    print OUT ">$subcon_name\n";
    print OUT $seq, "\n";
}
$seq = '';
close (OUT);
close (FILE);

my @target;
$seq = '';
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
        undef @sub_target;
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
        system ("nucmer --maxmatch -p $prefix_out -l $min_nucmer_match -banded -D 5 $ref_file $subcontig_file");
    }
    elsif (@target_subfiles > 1){
        foreach my $subfile (@target_subfiles){
            my $file_prefix = $1 if ($subfile =~ /(\S+)\.fasta$/);
            system ("nucmer --maxmatch -p $file_prefix -l $min_nucmer_match -banded -D 5 $subfile $subcontig_file");
            push @nucmer_subfiles, "$file_prefix.delta";
            print STDERR "Nucmer completed ($file_prefix.delta)\n";
        }
    }
}


if (@nucmer_subfiles > 0){
    open (OUT, "> $prefix_out.delta");
    print OUT "$ref_file $subcontig_file\n";
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

sub nucmer_align{
    my $subfile = shift;
    my $file_prefix = $1 if ($subfile =~ /(\S+)\.fasta$/);
    system ("nucmer --maxmatch -p $file_prefix -l $min_nucmer_match -banded -D 5 $subfile $subcontig_file 2> /dev/null");
    threads->yield();
    sleep 1;
    return ("$file_prefix.delta");
}


system ("delta-filter -o 95 -i $min_identity $prefix_out.delta > $prefix_out.fdelta");
system ("show-coords -lrcT $prefix_out.fdelta | sort -k13 -k1n -k2n > $prefix_out.coords");
#system ("show-tiling -l 1 -i 0 -V 0 $prefix_out.fdelta > $prefix_out.tilling");
if ($error_correct == 0){
    system ("$Bin/correct_scafs_coords.pl -q $scaf_file -s $subcontig_file -a $prefix_out.coords -p $prefix_out -mi $min_identity -ma $min_align_len -ml $min_subcon_len -mg $max_gap_len -is $max_indel_size -r $ref_file");
}
else{
    system ("$Bin/correct_scafs_coords.pl -q $scaf_file -s $subcontig_file -a $prefix_out.coords -p $prefix_out -mi $min_identity -ma $min_align_len -ml $min_subcon_len -mg $max_gap_len -is $max_indel_size -r $ref_file -e");
}
