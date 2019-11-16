#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

my $reference_file = '';         # -r
my $limited_mismatch_num = 2;	 # -n
my $limited_mismatch_rate = 0;	 # -m
my $out_prefix = 'coval-filter'; # -p
my $help;			 # -h
my $flag_type = 'BWA';           # -b   this is automatically specified

my $read_length = 0;


GetOptions(
    'ref=s' => \$reference_file,
    'num=i' => \$limited_mismatch_num,
    'mrate=f' => \$limited_mismatch_rate,
    'pref=s' => \$out_prefix,
    'bwa=s' => \$flag_type,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;
pod2usage(-verbose => 0) unless (-e $reference_file);
pod2usage(-verbose => 0) if (($limited_mismatch_num != 2) and ($limited_mismatch_rate != 0));

unless (@ARGV){
print "Your sam file name should be added to the argument.\n";
exit;
}
my $input_file = shift (@ARGV);
if ($input_file eq '-'){
}
elsif (!-f $input_file){
    die "inputfile does not exist: $!\n";
}


=head1 SYNOPSIS

 coval filter [options] <input_sorted_bam/sam_file> (Read names and reference sequence names in sam file should not contain characters '||' or '='.)
  Output:
   out_prefix.bam/sam

  Options:
   --ref or -r <STR>	   reference fasta file used for the alignment
   --pref or -p <STR>      prefix of output file
   --num or -n <INT>       maximum number of mismatches contained in a read [default: 2] (incompatible with --mrate)
   --mrate or -m <FLOAT>   maximum rate of mismatches contained in a read [0..1.0] (incompatible with --num)
   --help or -h            output help message

=cut


my %chr_seq;
open (FILE, $reference_file) or die "$!";
    my $ref_chr_name;
    my $count_chr = 0;
    my $chrseq = '';
    while (my $line = <FILE>){
        chomp $line;
	if ($line =~ /^>(\S*)/){
	    if ($count_chr > 0){
		$chr_seq{$ref_chr_name} = $chrseq;
		$chrseq = '';
	    }
	    $ref_chr_name = $1;
	    if (($ref_chr_name =~ /==/) or ($ref_chr_name =~ /\|\|/)){
		die "Reference chromosome/contig names contain unacceptable characters (== or \|\|): $!";
	    }
	    $count_chr ++;
	}
	else{
	    $chrseq .= uc $line;
	}
    }
    $chr_seq{$ref_chr_name} = $chrseq;
    $chrseq = '';
close (FILE);

my $terminal_mismatch = 0;
my $read_number = 0;
my $read_num_after_filt = 0;

my ($m1, $m2, $m3, $m4, $m5, $m6, $m7, $m8, $m9, $m10, $m10o) = (0,0,0,0,0,0,0,0,0,0,0);

open (FILE, $input_file) or die "$input_file: $!";
    while (my $line = <FILE>){
	if ($line =~ /^\@/){
	    print $line;
	    next;
	}
        chomp $line;
        my @line = split (/\s+/, $line);
	next if (@line < 11);
	next if ($line[5] =~ /\*/);
	next if ($line[1] == 4) or ($line[1] == 20) or ($line[1] == 7);
	$read_number++;
	if (($line[5] =~ /H/) or ($line[5] =~ /^\d+S.+S$/) or ($line[5] =~ /^\d+S\d+M\d+[DI]/) or ($line[5] =~ /\d+[DI]\d+M\d+S$/) or ($line[5] =~ /\d+[DI]\d+M\d+[DI]\d+M\d+[DI]/)){
	    next;
	}
        $line[2] =~ /^(\S*)/;
        my $chrom = $1;
        my $pos = $line[3];
	next if ($line[5] =~ /^\d+S/) and ($pos > 2);
        if ($line[5] =~ /^(\d+)M\d+S$/){
            if ($pos + $1 < length $chr_seq{$chrom}){
                next;
            }
        }
	$read_length = length $line[9];
	$line[9] = uc $line[9];
	my $read_name = $line[0];
	if ($read_name =~ /^[^:]+:([^#]+)#0/){
	    $read_name = $1;
	}
	elsif ($read_name =~ /^[^:]+:([^\/]+)\//){
	    $read_name = $1;
	}
	elsif ($read_name =~ /^[^:]+:(\S+)/){
	    $read_name = $1;
	}
	elsif ($read_name =~ /_*([^\/_]+)\//){
            $read_name = $1;
        }
	if ($limited_mismatch_rate != 0){
	    $limited_mismatch_num = int ($read_length * $limited_mismatch_rate + 0.5);
	}

	if ($line[5] eq $read_length . 'M'){
	    my $ref_seq = substr ($chr_seq{$chrom}, $line[3] - 1, $read_length);
	    if (length ($ref_seq) < $read_length){
		next;
	    }
	    if ($line[9] eq $ref_seq){
		print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		$read_num_after_filt ++;
	    }
	    else{
		my $total_mismatch = &count_mismatch_T ($line[9], $ref_seq);
		if ($total_mismatch <= $limited_mismatch_num){
		    print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		    $read_num_after_filt ++;
		}
	    }
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)I(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    next if (length ($chr_seq{$chrom}) - $line[3] - $1 + 1 < $3);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $3);
	    my $ref12 = $ref1 . $ref2;
	    next if (length $ref12 != $1 + $3);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1 + $2, $3);
	    my $read12 = $read1 . $read2;
	    if ($read12 eq $ref12){
		print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		$read_num_after_filt ++;
	    }
	    elsif ($limited_mismatch_num >= 2){
		my $total_mismatch = &count_mismatch_T ($read12, $ref12);
		if ($total_mismatch <= $limited_mismatch_num - 1){
		    print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		    $read_num_after_filt ++;
		}
	    }
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)D(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    next if (length ($chr_seq{$chrom}) - $line[3] - $1 -$2 + 1 < $3);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $3);
	    my $ref12 = $ref1 . $ref2;
	    next if (length $ref12 != $1 + $3);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1, $3);
	    my $read12 = $read1 . $read2;
	    if ($read12 eq $ref12){
		print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		$read_num_after_filt ++;
	    }
	    elsif ($limited_mismatch_num >= 2){
		my $total_mismatch = &count_mismatch_T ($read12, $ref12);
		if ($total_mismatch <= $limited_mismatch_num - 1){
		    print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		    $read_num_after_filt ++;
		}
	    }
	}
	elsif (($line[5] =~ /(\d+)M(\d+)([ID])(\d+)M(\d+)([ID])(\d+)M/) and ($limited_mismatch_num >= 2)){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $4) if ($3 eq 'I');
	    $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $4) if ($3 eq 'D');
	    
	    my $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $4 - 1, $7) if (($3 eq 'I') and ($6 eq 'I'));
	    $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $4 + $5 - 1, $7) if (($3 eq 'I') and ($6 eq 'D'));
	    $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $4 - 1, $7) if (($3 eq 'D') and ($6 eq 'I'));
	    $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $4 + $5 - 1, $7) if (($3 eq 'D') and ($6 eq 'D'));
	    next unless (defined $ref3);
	    my $ref12 = $ref1 . $ref2 . $ref3;
	    
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1 + $2, $4) if ($3 eq 'I');
	    $read2 = substr ($line[9], $1, $4) if ($3 eq 'D');
	    next if (length $ref12 != $1 + $4 + $7);
	    my $read3 = substr ($line[9], $1 + $2 + $4 + $5, $7) if (($3 eq 'I') and ($6 eq 'I'));
	    $read3 = substr ($line[9], $1 + $2 + $4, $7) if (($3 eq 'I') and ($6 eq 'D'));
	    $read3 = substr ($line[9], $1 + $4 + $5, $7) if (($3 eq 'D') and ($6 eq 'I'));
	    $read3 = substr ($line[9], $1 + $4, $7) if (($3 eq 'D') and ($6 eq 'D'));
	    my $read12 = $read1 . $read2 . $read3;	    
	    if ($read12 eq $ref12){
		print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		$read_num_after_filt ++;
	    }
	    elsif ($limited_mismatch_num >= 3){
		my $total_mismatch = &count_mismatch_T ($read12, $ref12);
		if ($total_mismatch <= $limited_mismatch_num - 2){
		    print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		    $read_num_after_filt ++;
		}
	    }
	}
	elsif ($line[5] =~ /^(\d+)S(\d+)M$/){
	    my $ref = substr ($chr_seq{$chrom}, $line[3] - 1, $2);
	    my $read = substr ($line[9], $1, $2);
	    if ($read eq $ref){
		print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		$read_num_after_filt ++;
	    }
	    elsif ($limited_mismatch_num >= 2){
		my $total_mismatch = &count_mismatch_T ($read, $ref);
		if ($total_mismatch <= $limited_mismatch_num - 1){
		    print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		    $read_num_after_filt ++;
		}
	    }
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)S$/){
	    my $ref = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $read = substr ($line[9], 0, $1);
	    if ($read eq $ref){
		print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		$read_num_after_filt ++;
	    }
	    elsif ($limited_mismatch_num >= 2){
		my $total_mismatch = &count_mismatch_T ($read, $ref);
		if ($total_mismatch <= $limited_mismatch_num - 1){
		    print "$read_name\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\n";
		    $read_num_after_filt ++;
		}
	    }
	}
    }
    print '@<<EOF>>', "\n";
close (FILE);

=pod
print STDERR "\n<< Number of mismatch reads before filtering >>\n";
print STDERR "1-mismatch reads = ", $m1, "\n";
print STDERR "2-mismatch reads = ", $m2, "\n";
print STDERR "3-mismatch reads = ", $m3, "\n";
print STDERR "4-mismatch reads = ", $m4, "\n";
print STDERR "5-mismatch reads = ", $m5, "\n";
print STDERR "6-mismatch reads = ", $m6, "\n";
print STDERR "7-mismatch reads = ", $m7, "\n";
print STDERR "8-mismatch reads = ", $m8, "\n";
print STDERR "9-mismatch reads = ", $m9, "\n";
print STDERR "10-mismatch reads = ", $m10, "\n";
print STDERR ">10-mismatch reads = ", $m10o, "\n";
=cut

print STDERR "Total number of aligned reads: $read_number\n";
print STDERR "Total number of reads after filtering: $read_num_after_filt\n";


###############################################################################

sub count_mismatch_T {
    my ($read, $ref) = @_;
    my $count = 0;
    return (100) if (length $read != length $ref);
    for (my $i = 0; $i < length $read; $i++){
	my $read_i = substr ($read, $i, 1);
	my $ref_i = substr ($ref, $i, 1);
	if ($read_i ne $ref_i){
	    $count++;
	}
    }
=pod
    $m1++ if ($count == 1);
    $m2++ if ($count == 2);
    $m3++ if ($count == 3);
    $m4++ if ($count == 4);
    $m5++ if ($count == 5);
    $m6++ if ($count == 6);
    $m7++ if ($count == 7);
    $m8++ if ($count == 8);
    $m9++ if ($count == 9);
    $m10++ if ($count == 10);
    $m10o++ if ($count > 10);
=cut
    return ($count);
}
