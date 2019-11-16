#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long qw(:config posix_default no_ignore_case);
use Pod::Usage;
use File::Basename;

my $target_scaf_file = '';     # -t
my $min_gap_size = 1;       # -mg
my $min_subcon_size = 150;  # -ms
my $out_prefix = '';        # -p

GetOptions(
    'target_scaf|t=s' => \$target_scaf_file,
    'min_subcon|ms=i' => \$min_subcon_size,
    'min_gap_size|mg=i' => \$min_gap_size,
    'prefix|p=s' => \$out_prefix,
) or pod2usage(-verbose => 0);

my %target_subcon_seq;
my $target_name = '-';
my $Tseq = '';

open (FILE, $target_scaf_file) or die "$target_scaf_file is not found: $!\n";
while (my $line = <FILE>){
    $line =~ s/(\r\n|\n|\r)//g;
    if ($line =~ /^>/){
        if ($target_name ne '-'){
            $Tseq =~ s/^[Nn]+// if ($Tseq =~ /^[Nn]/);
            $Tseq =~ s/[Nn]+$// if ($Tseq =~ /[Nn]$/);
            if ($Tseq =~ /[Nn]{$min_gap_size}/){
                &assign_Tsubcont ($target_name, $Tseq);
            }
            else{
                my $count_subcon_04d = sprintf ("%04d", 1);
                my $subcon_name = $target_name . '/' . $count_subcon_04d;
                $target_subcon_seq{$subcon_name} = $Tseq;
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
if ($Tseq =~ /[Nn]{$min_gap_size}/){
    &assign_Tsubcont ($target_name, $Tseq);
}
else{
    my $count_subcon_04d = sprintf ("%04d", 1);
    my $subcon_name = $target_name . '/' . $count_subcon_04d;
    $target_subcon_seq{$subcon_name} = $Tseq;
}
$Tseq = '';
close (FILE);

my $target_subcon_file = "$out_prefix-subcon.fa";

open (OUT, "> $target_subcon_file");
foreach my $name (sort keys %target_subcon_seq){
    print OUT ">$name\n";
    print OUT $target_subcon_seq{$name}, "\n";
}
close (OUT);

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
        $carry_gap_len = 0;
    }
    if ($seq =~ /$subseq([A-Z]+?)$/){
        $count_subcon ++;
        my $count_subcon_04d = sprintf ("%04d", $count_subcon);
        my $subcon_name = $name . '/' . $count_subcon_04d;
        $target_subcon_seq{$subcon_name} = $1;
    }
}