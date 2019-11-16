#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

# connect subcontigs from gmcloser-blast-LR multi-thread jobs

my $subfile = '';
my $query_file = '';
my $num = 0;
my $out_prefix = '';
my $temp_dir = '';
my $gapA_file = '';
my $gap5_file = '';
my $gap3_file = '';
my $pe_align_tag = 2;
my $max_insert_size = 0;
my $read_length = 0;
my $T_PEalignF_file = '';
my $T_PEalignR_file = '';
my $Q_PEalignF_file = '';
my $Q_PEalignR_file = '';
my $min_identity = 99;
my $min_match_length = 300;
my $min_match_len_local = 20;
my $connect = 0;
my $ite_num = 1;


GetOptions(
    'subfile|s=s' => \$subfile,
    'query|q=s' => \$query_file,
    'num|n=i' => \$num,
    'prefix|p=s' => \$out_prefix,
    'tmpdir|td=s' => \$temp_dir,
    'gapa|ga=s' => \$gapA_file,
    'gapu|gu=s' => \$gap5_file,
    'gapd|gd=s' => \$gap3_file,
    'tag|t=i' => \$pe_align_tag,
    'max_insert|i=i' => \$max_insert_size,
    'min_ident|mi=i' => \$min_identity,
    'min_match_len|mm=i' => \$min_match_length,
    'min_len|ml=i' => \$min_match_len_local,
    'readlen|r=i' => \$read_length,
    'tf_file|tf=s' => \$T_PEalignF_file,
    'tr_file|tr=s' => \$T_PEalignR_file,
    'qf_file|qf=s' => \$Q_PEalignF_file,
    'qr_file|qr=s' => \$Q_PEalignR_file,
    'connect|cn' => \$connect,
    'ite_num|in=i' => \$ite_num
);

die "subfile is not specified: $!\n" if ($subfile eq '');
die "query_file is not specified: $!\n" if ($query_file eq '');
die "out_prefix is not specified: $!\n" if ($out_prefix eq '');
die "temp_dir is not specified: $!\n" if ($temp_dir eq '');
die "gapA_file is not specified: $!\n" if ($gapA_file eq '');
die "gap5_file is not specified: $!\n" if ($gap5_file eq '');
die "gap3_file is not specified: $!\n" if ($gap3_file eq '');
die "PE_align_tag should be 0 or 1: $!\n" if ($pe_align_tag != 0) and ($pe_align_tag != 1);
die "max_insert_size is not specified: $!\n" if ($max_insert_size == 0) and ($pe_align_tag == 1);
die "read_length is not specified: $!\n" if ($read_length == 0) and ($pe_align_tag == 1);
die "T_PEalignF_file is not specified: $!\n" if ($T_PEalignF_file eq '') and ($pe_align_tag == 1);
die "T_PEalignR_file is not specified: $!\n" if ($T_PEalignR_file eq '') and ($pe_align_tag == 1);
die "Q_PEalignF_file is not specified: $!\n" if ($Q_PEalignF_file eq '') and ($pe_align_tag == 1);
die "Q_PEalignR_file is not specified: $!\n" if ($Q_PEalignR_file eq '') and ($pe_align_tag == 1);

my $connect_file = "$temp_dir/$out_prefix-subcon-connect-$num.txt";
my $delete_tqname_file = "$temp_dir/$out_prefix-delete_tqname-$num.txt";
my $targetname = '';
my $gaplen = 0;
my $seq = '';
my $connect_num_local = 0;
my %target_subcon_seq_2;
my %target_subcon_num;
my %target_gap_size;
my %query_seq;
my %gap_A;
my %gap_5;
my %gap_3;
my %T_PEalign_F;
my %T_PEalign_R;
my %Q_PEalign_F;
my %Q_PEalign_R;

open (FILE, $subfile) or die "$subfile is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^>/){
        if ($seq ne ''){
            my $tname_base = $1 if ($targetname =~ /(\S+)\/(\d+)$/);
            my $sub_no = $2;
            $sub_no =~ s/^0*//;
            ${$target_subcon_seq_2{$tname_base}}{$sub_no} = $seq;
            $target_gap_size{$targetname} = $gaplen;
            $seq = '';
        }
        ($targetname, $gaplen) = split (/\s+/, $line);
        $targetname = $1 if ($targetname =~ /^>(.+)/);
        $gaplen = 0 if (!defined $gaplen) or ($gaplen eq '');
    }
    else{
        $seq .= $line;
    }
}
my $tname_base = $1 if ($targetname =~ /(\S+)\/(\d+)$/);
my $sub_no = $2;
$sub_no =~ s/^0*//;
${$target_subcon_seq_2{$tname_base}}{$sub_no} = $seq;
$target_gap_size{$targetname} = $gaplen;
$seq = '';
close (FILE);

foreach my $tbase (keys %target_subcon_seq_2){
    foreach my $sub_no (sort {$b <=> $a} keys %{$target_subcon_seq_2{$tbase}}){
        $target_subcon_num{$tbase} = $sub_no;
        last;
    }
}

my $query_name = '';
my $Qseq = '';
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

if ($pe_align_tag == 1){
    open (FILE, $T_PEalignF_file) or die "$T_PEalignF_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($tname, $pos, $reads) = split (/\t/, $line);
        ${$T_PEalign_F{$tname}}{$pos} = $reads;
    }
    close (FILE);
    open (FILE, $T_PEalignR_file) or die "$T_PEalignR_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($tname, $pos, $reads) = split (/\t/, $line);
        ${$T_PEalign_R{$tname}}{$pos} = $reads;
    }
    close (FILE);
    open (FILE, $Q_PEalignF_file) or die "$Q_PEalignF_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($qname, $pos, $reads) = split (/\t/, $line);
        ${$Q_PEalign_F{$qname}}{$pos} = $reads;
    }
    close (FILE);
    open (FILE, $Q_PEalignR_file) or die "$Q_PEalignR_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($qname, $pos, $reads) = split (/\t/, $line);
        ${$Q_PEalign_R{$qname}}{$pos} = $reads;
    }
    close (FILE);
}

open (FILE, $gapA_file) or die "$gapA_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($tname, $gap_info) = split (/\t/, $line);
    $gap_A{$tname} = $gap_info;
}
close (FILE);

open (FILE, $gap5_file) or die "$gap5_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($tname, $gap_info) = split (/\t/, $line);
    $gap_5{$tname} = $gap_info;
}
close (FILE);

open (FILE, $gap3_file) or die "$gap3_file is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($tname, $gap_info) = split (/\t/, $line);
    $gap_3{$tname} = $gap_info;
}
close (FILE);


open (OUT1, "> $connect_file");
open (OUT2, "> $delete_tqname_file");
foreach my $tbase (keys %target_subcon_seq_2){
    foreach my $subno (sort {$a <=> $b} keys %{$target_subcon_seq_2{$tbase}}){
        my $tname = $tbase. '/' . sprintf ("%04d", $subno);
        my $down_tname = $tbase. '/' . sprintf ("%04d", $subno + 1);
        next if ($subno == $target_subcon_num{$tbase});       
        my $gap_5_unmatch_len = 0;
        my $gap_3_unmatch_len = 0;
        my $gap5_seq = '';
        my $gap3_seq = '';
        my $gap5_qname = '';
        my $gap3_qname = '';
        my $gap5_TQtag = '';
        my $gap3_TQtag = '';
        if (!exists $gap_A{$tname}){
            if ((exists $gap_5{$tname}) or (exists $gap_3{$down_tname})){
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
                    $gap5_seq = ${$target_subcon_seq_2{$tbase}}{$subno};
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
                    $gap3_seq = ${$target_subcon_seq_2{$tbase}}{$subno + 1};
                    $gap3_TQtag = 'T5';
                }
                my ($match, $match_len1, $match_len2, $identity2) = &local_align (\$gap5_seq, \$gap3_seq);
                my $rate_PE = 0;
                my $read_num = 0;
                if (($match == 1) and ($pe_align_tag == 1)){
                    my $gap3_len = length $gap3_seq;
                    my $gap5_len = length $gap5_seq;
                    next if ($gap5_len < $read_length) or ($gap3_len < $read_length);
                    my $gap5_name1 = $tname if ($gap5_TQtag eq 'T3');
                    $gap5_name1 = $gap5_qname if ($gap5_TQtag eq 'Q5') or ($gap5_TQtag eq 'Q3');
                    my $gap3_name1 = $down_tname if ($gap3_TQtag eq 'T5');
                    $gap3_name1 = $gap3_qname if ($gap3_TQtag eq 'Q5') or ($gap3_TQtag eq 'Q3');
                    ($rate_PE, $read_num) = &find_aligned_read_connect ($gap5_name1, $gap5_TQtag, $gap5_len, $gap3_name1, $gap3_TQtag, $gap3_len, $match_len1);
                }
                if ((($match == 1) and ($pe_align_tag == 1) and (($match_len1 >= $min_match_length) or (($match_len1 >= $min_match_len_local) and ($identity2 >= $min_identity) and ($read_num < 10)) or ($rate_PE >= 0.08))) or 
                (($match == 1) and ($pe_align_tag == 0) and ($match_len1 >= $min_match_len_local) and ($identity2 >= $min_identity))){
                    if (($gap5_qname ne '') and ($gap3_qname ne '')){
                        if ($match_len1 >= $gap_5_unmatch_len + $gap_3_unmatch_len){
                            my $gap_size = $target_gap_size{$tname};
                            next if ($gap_size >= 500) and ($rate_PE < 0.08);
                            my $overlen = $match_len1 - $gap_5_unmatch_len - $gap_3_unmatch_len;
                            print OUT1 "$tname\tC\t$overlen\n";
                        }
                        elsif (($match_len1 <= $gap_5_unmatch_len) and ($match_len2 <= $gap_3_unmatch_len)){
                            my $gap_seq_1 = substr ($gap5_seq, -$gap_5_unmatch_len, $gap_5_unmatch_len);
                            substr ($gap_seq_1, -$match_len1, $match_len1, '');
                            my $gap_seq_2 = substr ($gap3_seq, 0, $gap_3_unmatch_len);
                            my $gap_seq = $gap_seq_1 . $gap_seq_2;
                            print OUT1 "$tname\tA\t$gap_seq\t$gap5_qname\t$gap3_qname\n";
                        }
                        elsif ($match_len1 <= $gap_5_unmatch_len){
                            my $gap_seq = substr ($gap5_seq, -$gap_5_unmatch_len, $gap_5_unmatch_len);
                            substr ($gap_seq, -($match_len2 - $gap_3_unmatch_len), $match_len2 - $gap_3_unmatch_len, '');
                            print OUT1 "$tname\tA\t$gap_seq\t$gap5_qname\t$gap3_qname\n";
                        }
                        elsif ($match_len2 <= $gap_3_unmatch_len){
                            my $gap_seq = substr ($gap3_seq, 0, $gap_3_unmatch_len);
                            substr ($gap_seq, 0, $match_len1 - $gap_5_unmatch_len, '');
                            print OUT1 "$tname\tA\t$gap_seq\t$gap5_qname\t$gap3_qname\n";
                        }
                    }
                    elsif ($gap5_qname ne ''){
                        if ($match_len1 >= $gap_5_unmatch_len){
                            my $gap_size = $target_gap_size{$tname};
                            next if ($gap_size >= 500) and ($rate_PE < 0.08);
                            my $overlen = $match_len1 - $gap_5_unmatch_len;
                            print OUT1 "$tname\tC\t$overlen\n";
                        }
                        else{
                            my $gap_seq = substr ($gap5_seq, -$gap_5_unmatch_len, $gap_5_unmatch_len);
                            substr ($gap_seq, -$match_len1, $match_len1, '');
                            print OUT1 "$tname\tA\t$gap_seq\t$gap5_qname\n";
                        }
                    }
                    elsif ($gap3_qname ne ''){
                        if ($match_len2 >= $gap_3_unmatch_len){
                            my $gap_size = $target_gap_size{$tname};
                            next if ($gap_size >= 500) and ($rate_PE < 0.08);
                            my $overlen = $match_len2 - $gap_3_unmatch_len;
                            print OUT1 "$tname\tC\t$overlen\n";
                        }
                        else{
                            my $gap_seq = substr ($gap3_seq, 0, $gap_3_unmatch_len);
                            substr ($gap_seq, 0, $match_len2, '');
                            print OUT1 "$tname\tA\t$gap_seq\t$gap3_qname\n";
                        }
                    }
                }
                elsif ($match == 0){
                    my $gap_size = $target_gap_size{$tname};
                    if ($gap_5_unmatch_len + $gap_3_unmatch_len > $gap_size * 1.5){
                        if ($gap5_qname ne ''){
                            print OUT2 "$tname\t$gap5_qname\n";
                        }
                        if ($gap3_qname ne ''){
                            print OUT2 "$tname\t$gap3_qname\n";
                        }
                    }
                }
            }
            else{
                if (($connect == 1) and ($ite_num == 1)){
                    $gap5_seq = ${$target_subcon_seq_2{$tbase}}{$subno};
                    $gap5_TQtag = 'T3';
                    $gap3_seq = ${$target_subcon_seq_2{$tbase}}{$subno + 1};
                    $gap3_TQtag = 'T5';
                    next if (length $gap5_seq < 20) or (length $gap3_seq < 20);
                    my ($match, $match_len1, $match_len2, $identity2) = &local_align (\$gap5_seq, \$gap3_seq);
                    my $rate_PE = 0;
                    my $read_num = 0;
                    if (($match == 1) and ($pe_align_tag == 1)){
                        my $gap3_len = length $gap3_seq;
                        my $gap5_len = length $gap5_seq;
                        ($rate_PE, $read_num) = &find_aligned_read_connect ($tname, $gap5_TQtag, $gap5_len, $down_tname, $gap3_TQtag, $gap3_len, $match_len1);
                    }
                    my $gap_size = $target_gap_size{$tname};
                    next if ($gap_size >= 500) and ($rate_PE < 0.08);
                    next if ($identity2 < 98) and ($rate_PE < 0.08);
                    if (($match_len1 >= $min_match_length) or (($pe_align_tag == 1) and ($match_len1 >= $min_match_len_local) and ($identity2 >= $min_identity) and ($read_num < 10)) or ($rate_PE >= 0.08)){
                        print OUT1 "$tname\tC\t$match_len1\n";
                    }
                }
            }
        }
    }
}
close (OUT1);
close (OUT2);

sub local_align{                    # align between two query sequences
    my ($ref_seq1, $ref_seq2) = @_;        # $qseq1: upstream subcon + unamtch region of the overlapped query, $qseq2: unmatch region of query that overlaps target subcon
    my $seq1 = ${$ref_seq1};
    my $seq2 = ${$ref_seq2};
    open (NEWFILE1, "> temp/$out_prefix-$num-seq1.fa");
    print NEWFILE1 '>seq1', "\n", $seq1;
    close (NEWFILE1);
                
    open (NEWFILE1, "> temp/$out_prefix-$num-seq2.fa");
    print NEWFILE1 '>seq2', "\n", $seq2;
    close (NEWFILE1);
    my $mutation_rate = 3;
    my @result = `yass -O 5 -m $mutation_rate -r 0 temp/$out_prefix-$num-seq1.fa temp/$out_prefix-$num-seq2.fa 2>/dev/null`;
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
                if (exists ${$T_PEalign_F{$up_cont}}{$i}){
                    push @up_read_F1, split (/=/, ${$T_PEalign_F{$up_cont}}{$i});
                }
            }
            elsif ($up_tag eq 'Q3'){
                if (exists ${$Q_PEalign_F{$up_cont}}{$i}){
                    push @up_read_F1, split (/=/, ${$Q_PEalign_F{$up_cont}}{$i});
                }
            }
            elsif ($up_tag eq 'Q5'){
                if (exists ${$Q_PEalign_R{$up_cont}}{$i}){
                    push @up_read_R1, split (/=/, ${$Q_PEalign_R{$up_cont}}{$i});
                }
            }
        }
        for (my $i = $scan_start1; $i <= $scan_end1; $i++){
            if ($tag eq 'T5'){
                if (exists ${$T_PEalign_R{$cont}}{$i}){
                    push @read_R1, split (/=/, ${$T_PEalign_R{$cont}}{$i});
                }
            }
            elsif ($tag eq 'Q5'){
                if (exists ${$Q_PEalign_R{$cont}}{$i}){
                    push @read_R1, split (/=/, ${$Q_PEalign_R{$cont}}{$i});
                }
            }
            elsif ($tag eq 'Q3'){
                if (exists ${$Q_PEalign_F{$cont}}{$i}){
                    push @read_F1, split (/=/, ${$Q_PEalign_F{$cont}}{$i});
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
                if (exists ${$T_PEalign_F{$up_cont}}{$i}){
                    push @up_read_F2, split (/=/, ${$T_PEalign_F{$up_cont}}{$i});
                }
            }
            elsif ($up_tag eq 'Q3'){
                if (exists ${$Q_PEalign_F{$up_cont}}{$i}){
                    push @up_read_F2, split (/=/, ${$Q_PEalign_F{$up_cont}}{$i});
                }
            }
            elsif ($up_tag eq 'Q5'){
                if (exists ${$Q_PEalign_R{$up_cont}}{$i}){
                    push @up_read_R2, split (/=/, ${$Q_PEalign_R{$up_cont}}{$i});
                }
            }
        }
        for (my $i = $scan_start2; $i <= $scan_end2; $i++){
            if ($tag eq 'T5'){
                if (exists ${$T_PEalign_R{$cont}}{$i}){
                    push @read_R2, split (/=/, ${$T_PEalign_R{$cont}}{$i});
                }
            }
            elsif ($tag eq 'Q5'){
                if (exists ${$Q_PEalign_R{$cont}}{$i}){
                    push @read_R2, split (/=/, ${$Q_PEalign_R{$cont}}{$i});
                }
            }
            elsif ($tag eq 'Q3'){
                if (exists ${$Q_PEalign_F{$cont}}{$i}){
                    push @read_F2, split (/=/, ${$Q_PEalign_F{$cont}}{$i});
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
