#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

# Generate a contig fasta file whose contigs are broken at misassembled positions
# The misassembled positions are obtained from 'out.coords' file

my $contig_file = '';
my $align_pos_file = '';
my $ref_file = '';
my $min_identity = 97;
my $min_align_len = 200;
my $min_contig_len = 200;
my $min_coverage = 99;
my $prefix_out = 'out';
my $error_correct = 0;
my $max_indel_size = 100;
my $help;

GetOptions(
    'query|q=s' => \$contig_file,
    'align|a=s' => \$align_pos_file,
    'ref|r=s' => \$ref_file,
    'min_id|mi=i' => \$min_identity,
    'min_len|ml=i' => \$min_contig_len,
    'min_align|ma=i' => \$min_align_len,
    'min_cov|mc=i' => \$min_coverage,
    'prefix|p=s' => \$prefix_out,
    'error_correct|e' => \$error_correct,
    'max_indel|is=i' => \$max_indel_size,
    'help' => \$help
) or pod2usage(-verbose => 0);

pod2usage(-verbose => 0) if $help;
pod2usage(-verbose => 0) if ($contig_file eq '');
pod2usage(-verbose => 0) if ($align_pos_file eq '');
pod2usage(-verbose => 0) if ($ref_file eq '') and ($error_correct == 1);

=head1 SYNOPSIS

  gmvalue ver. 1.3
  
  Options:
   --query or -q <STR>      input contig fasta file (e.g., contig1.fa)
   --ref or -r <STR>        input reference file (e.g., ref.fa)
   --align or -a <STR>      input coords file from Nucmer outputs (e.g., align.coords)
   --min_id or -mi <INT>    minimum alignment identity (%) [default: 97]
   --min_cov or -mc <INT>   minimum coverage (%) of query (contig) aligned to a reference [default: 99]
   --min_align or -ma <INT> minimum alignment overlap length with the maximum allowable size of indels [default: 200]
   --min_len or -ml <INT>   minimum contig length to be considered [default: 200]
   --prefix or -p <STR>     prefix name of outputs
   --error_correct or -e    output an error-corrected contig set [default: false]
   --max_indel or -is <INT> maximum allowable size of indels (or distance between break points) [default: 100]
   --thread or -n           number of threads to run [default: 1]
   --help or -h             output help message
   
=cut

my %ref_seq;
my $chr_name = '-';
my $Rseq = '';
open (FILE, $ref_file) or die "$!";
while (my $line = <FILE>){
    $line =~ s/(\r\n|\n|\r)//g;
    if ($line =~ /^>(\S+)/){
        if ($chr_name ne '-'){
            $ref_seq{$chr_name} = uc $Rseq;
            $Rseq = '';
        }
        $chr_name = $1;
    }
    else{
        $Rseq .= $line;
    }
}
$ref_seq{$chr_name} = uc $Rseq;
$Rseq = '';
close (FILE);

my $query_name = '-';
my $Qseq = '';
my %query_seq;
my %query_length;
open (FILE, $contig_file) or die "$!";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^>(\S+)/){
        if ($query_name ne '-'){
            $query_seq{$query_name} = uc $Qseq;
            $query_length{$query_name} = length $Qseq;
            $Qseq = '';
        }
        $query_name = $1;
    }
    else{
        $Qseq .= $line;
    }
}
$query_seq{$query_name} = uc $Qseq;
$query_length{$query_name} = length $Qseq;
$Qseq = '';
close (FILE);

my %partial_match;
my %broken_query_seq;
my %broken_target_seq;
my %misassemble_query;
my %assemble_query;

my $middle_indel = 0;
my $local_error = 0;
my $translocate = 0;
my $inversion = 0;
my $small_cov = 0;
my $local_error_L = 0;

open (FILE, $align_pos_file) or die "$!";
while (my $line = <FILE>){
    chomp $line;
    next if ($line !~ /^\d+/);
    my @line = split (/\s+/, $line);
    my $qname = $line[12];
    my $strand = 'Plus' if ($line[2] < $line[3]);
    $strand = 'Minus' if ($line[2] > $line[3]);
    my $qstart = $line[2] if ($strand eq 'Plus');
    $qstart = $line[3] if ($strand eq 'Minus');
    my $qend = $line[2] if ($strand eq 'Minus');
    $qend = $line[3] if ($strand eq 'Plus');
    if ($line[10] >= $min_coverage){
        $assemble_query{$qname} = "$qstart=$qend";
        next;
    }
    my $tstart = $line[0];
    my $tend = $line[1];
    my $identity = $line[6];
    my $overlen = $line[4];
    my $qlen = $line[8];
    my $tlen = $line[7];
    my $tname = $line[11];
    if (!exists $assemble_query{$qname}){
        ${${${$partial_match{$tname}}{$qname}}{$strand}}{$tstart} = "$tend==$qstart==$qend==$overlen==$identity";
    }
}
close (FILE);

foreach my $tname (keys %partial_match){
    my $tlen = length $ref_seq{$tname};
    foreach my $qname (keys %{$partial_match{$tname}}){
        next if (exists $assemble_query{$qname});
        my $qlen = length $query_seq{$qname};
        foreach my $strand (keys %{${$partial_match{$tname}}{$qname}}){
            my $pre_tend = 0;
            my $pre_tstart = 0;
            my $pre_qend = 0;
            my $pre_qstart = 0;
            my $count_align = 0;
            my $align_num = scalar keys %{${${$partial_match{$tname}}{$qname}}{$strand}};
            foreach my $tstart (sort {$a<=>$b} keys %{${${$partial_match{$tname}}{$qname}}{$strand}}){
                $count_align ++;
                my ($tend, $qstart, $qend, $overlen, $ident) = split (/==/, ${${${$partial_match{$tname}}{$qname}}{$strand}}{$tstart});
                if ($count_align == 1){                     # considerer the gap regions and the 5'-terminus of the reference
                    if ($strand eq 'Plus'){
                        if ($qstart > 10){
                            if ($tstart <= 2){
                                $qstart = 1;
                            }
                            elsif ($tstart > 10){
                                my $upstream_seq = substr ($ref_seq{$tname}, $tstart - $qstart, $qstart - 6) if ($tstart - $qstart >= 0);
                                $upstream_seq = substr ($ref_seq{$tname}, 0, $tstart - 6) if ($tstart - $qstart < 0);
                                if ($upstream_seq =~ /^N+$/){
                                    $qstart = 1;
                                }
                            }
                        }
                    }
                    else{
                        if ($qend <= $qlen - 10){
                            if ($tstart <= 2){
                                $qend = $qlen;
                            }
                            elsif ($tstart > 10){
                                my $upstream_seq = substr ($ref_seq{$tname}, $tstart - ($qlen - $qend + 1), $qlen - $qend - 5) if ($tstart - ($qlen - $qend + 1) >= 0);
                                $upstream_seq = substr ($ref_seq{$tname}, 0, $qlen - $qend - 5) if ($tstart - ($qlen - $qend + 1) < 0);
                                if ($upstream_seq =~ /^N+$/){
                                    $qend = $qlen;
                                }
                            }
                        }
                    }
                }
                if ($count_align == $align_num){            # considerer the gap regions and the 3'-terminus of the reference
                    if ($strand eq 'Plus'){
                        if ($qend <= $qlen - 10){
                            if ($tend >= $tlen - 2){
                                $qend = $qlen;
                            }
                            elsif ($tend <= $tlen - 10){
                                my $downstream_seq = substr ($ref_seq{$tname}, $tend + 4, $qlen - $qend - 5) if ($tend + $qlen - $qend <= $tlen);
                                $downstream_seq = substr ($ref_seq{$tname}, $tend + 4) if ($tend + $qlen - $qend > $tlen);
                                if ($downstream_seq =~ /^N+$/){
                                    $qend = $qlen;
                                }
                            }
                        }
                    }
                    else{
                        if ($qstart > 10){
                            if ($tend >= $tlen - 2){
                                $qstart = 1;
                            }
                            elsif ($tend <= $tlen - 10){
                                my $downstream_seq = substr ($ref_seq{$tname}, $tend + 4, $qstart - 5) if ($tend + $qstart <= $tlen);
                                $downstream_seq = substr ($ref_seq{$tname}, $tend + 4) if ($tend + $qstart > $tlen);
                                if ($downstream_seq =~ /^N+$/){
                                    $qstart = 1;
                                }
                            }
                        }
                    }
                }
                if ($pre_tend > 0){
                    if ($tend <= $pre_tend){
                        next;
                    }
                    if (($qend > $pre_qend) and ($qstart < $pre_qstart)){
                        $pre_tend = $tend;
                        $pre_tstart = $tstart;
                        $pre_qend = $qend;
                        $pre_qstart = $qstart;
                        next;
                    }
                    my $t_distance = abs ($tstart - $pre_tend);
                    my $q_distance = abs ($qstart - $pre_qend) if ($strand eq 'Plus');
                    $q_distance = abs ($qend - $pre_qstart) if ($strand eq 'Minus');
                    if (($q_distance <= $max_indel_size) and ($t_distance <= $max_indel_size)){
#                    if ($t_distance <= $max_indel_size){
                        $middle_indel ++;
                        $pre_tend = $tend;
                        $pre_qend = $qend if ($strand eq 'Plus');
                        $pre_qstart = $qstart if ($strand eq 'Minus');
                        next;
                    }
                    else{
                        my $total_overlen = $pre_qend - $pre_qstart + 1;
                        if ($total_overlen / $qlen * 100 >= $min_coverage){
                            $assemble_query{$qname} = "$pre_qstart=$pre_qend";
                        }
                        else{
                            if (!exists $assemble_query{$qname}){
                                if ($total_overlen >= $min_align_len){
                                    ${${${$broken_query_seq{$qname}}{$tname}}{$strand}}{$pre_qstart} = $total_overlen;
                                    ${${${$broken_target_seq{$qname}}{$tname}}{$strand}}{$pre_qstart} = "$pre_tstart=$pre_tend";
                                }
                                else{
                                    ${$misassemble_query{$qname}}{$pre_qstart} = $total_overlen;
                                }
                            }
                        }
                    }
                }
                $pre_tend = $tend;
                $pre_tstart = $tstart;
                $pre_qend = $qend;
                $pre_qstart = $qstart;
            }
            my $total_overlen = $pre_qend - $pre_qstart + 1;
            if ($total_overlen / $qlen * 100 >= $min_coverage){
                $assemble_query{$qname} = "$pre_qstart=$pre_qend";
            }
            else{
                if (!exists $assemble_query{$qname}){
                    if ($total_overlen >= $min_align_len){
                        ${${${$broken_query_seq{$qname}}{$tname}}{$strand}}{$pre_qstart} = $total_overlen;
                        ${${${$broken_target_seq{$qname}}{$tname}}{$strand}}{$pre_qstart} = "$pre_tstart=$pre_tend";
                    }
                    else{
                        ${$misassemble_query{$qname}}{$pre_qstart} = $total_overlen;
                    }
                }
            }
        }
    }
}

my %broken_query_seq_2;
my %broken_query_cov;
my %broken_query_uncov;
foreach my $qname (keys %broken_query_seq){
    next if (exists $assemble_query{$qname});
    foreach my $tname (keys %{$broken_query_seq{$qname}}){
        foreach my $strand (keys %{${$broken_query_seq{$qname}}{$tname}}){
            my %selected;
            my $sub_num = 0;
            foreach my $qstart (sort {${${${$broken_query_seq{$qname}}{$tname}}{$strand}}{$b} <=> ${${${$broken_query_seq{$qname}}{$tname}}{$strand}}{$a}} keys %{${${$broken_query_seq{$qname}}{$tname}}{$strand}}){
                my $qend = $qstart + ${${${$broken_query_seq{$qname}}{$tname}}{$strand}}{$qstart} - 1;
                my $overlap_flag = 0;
                if (scalar keys %selected >= 1){
                    foreach my $start2 (keys %selected){
                        my $end2 = $selected{$start2};
                        if ($qstart == $start2){
                            $overlap_flag = 1;
                            last;
                        }
                        elsif (($qstart >= $start2) and ($qend <= $end2)){
                            $overlap_flag = 1;
                            last;
                        }
                        elsif (($qstart < $start2) and ($qend > $start2) and (($qend - $start2 >= ($end2 - $start2) * 0.8) or ($qend - $start2 >= ($qend - $qstart) * 0.8))){
                            $overlap_flag = 1;
                            last;
                        }
                        elsif (($qstart < $end2) and ($qend > $end2) and (($end2 - $qstart >= ($end2 - $start2) * 0.8) or ($end2 - $qstart >= ($qend - $qstart) * 0.8))){
                            $overlap_flag = 1;
                            last;
                        }
        #                if ((($qstart > $start2) and ($qstart <= $end2 - $max_indel_size)) or (($qstart < $start2) and ($qend >= $start2 + $max_indel_size))){
        #                    $overlap_flag = 1;
        #                    last;
        #                }
                    }
                    if ($overlap_flag == 0){
                        $selected{$qstart} = $qend;
                    }
                }
                else{
                    $selected{$qstart} = $qend;
                }
            }
            my $qlen = length $query_seq{$qname};
            foreach my $qstart (keys %selected){
                ${${${$broken_query_seq_2{$qname}}{$tname}}{$strand}}{$qstart} = $selected{$qstart};
            }
            next if (!exists ${${$broken_query_seq_2{$qname}}{$tname}}{$strand});
            my $pre_qend = 0;
            my $overlap = 0;
            foreach my $qstart (sort {$a <=> $b} keys %{${${$broken_query_seq{$qname}}{$tname}}{$strand}}){
                my $qend = $qstart + ${${${$broken_query_seq{$qname}}{$tname}}{$strand}}{$qstart} - 1;
                if ($pre_qend == 0){
                    $overlap += $qend - $qstart + 1;
                }
                else{
                    if ($qstart <= $pre_qend){
                        $overlap += $qend - $pre_qend + 1;
                    }
                    else{
                        $overlap += $qend - $qstart + 1;
                    }
                }
                $pre_qend = $qend;
            }
            if ($overlap / $qlen >= $min_coverage * 0.01){
                if (!exists $broken_query_cov{$qname}){
                    $broken_query_cov{$qname} = "$tname==$strand==$overlap";
                }
                else{
                    my ($tname2, $strand2, $overlap2) = split (/==/, $broken_query_cov{$qname});
                    if ($overlap > $overlap2){
                        $broken_query_cov{$qname} = "$tname==$strand==$overlap";
                    }
                }
            }
            else{
                if (!exists $broken_query_uncov{$qname}){
                    $broken_query_uncov{$qname} = "$tname==$strand==$overlap";
                }
                else{
                    my ($tname2, $strand2, $overlap2) = split (/==/, $broken_query_uncov{$qname});
                    if ($overlap > $overlap2){
                        $broken_query_uncov{$qname} = "$tname==$strand==$overlap";
                    }
                }
            }
        }
    }
}

my %trans_pos;
my %inv_pos;
foreach my $qname (keys %broken_query_seq_2){
    my $qlen = length $query_seq{$qname};
    if (exists $broken_query_cov{$qname}){
        my ($tname2, $strand2, $overlap2) = split (/==/, $broken_query_cov{$qname});
        foreach my $tname (keys %{$broken_query_seq_2{$qname}}){
            if ($tname ne $tname2){
                delete ${$broken_query_seq_2{$qname}}{$tname};
                next;
            }
            foreach my $strand (keys %{${$broken_query_seq_2{$qname}}{$tname2}}){
                if ($strand ne $strand2){
                    delete ${${$broken_query_seq_2{$qname}}{$tname}}{$strand};
                    next;
                }
            }
        }
    }
    elsif (exists $broken_query_uncov{$qname}){
        my ($tname2, $strand2, $overlap2) = split (/==/, $broken_query_uncov{$qname});
        my @start_end;
        if ((scalar keys %{$broken_query_seq_2{$qname}} > 1) or (scalar keys %{${$broken_query_seq_2{$qname}}{$tname2}} > 1)){
            foreach my $qstart (sort {$a <=> $b} keys %{${${$broken_query_seq_2{$qname}}{$tname2}}{$strand2}}){
                my $qend = ${${${$broken_query_seq_2{$qname}}{$tname2}}{$strand2}}{$qstart};
                push @start_end, "$qstart=$qend";
            }
        }
        if (scalar keys %{$broken_query_seq_2{$qname}} > 1){
            my $pre_qend2 = 0;
            my $pre_qstart2 = 0;
            foreach my $range (@start_end){
                my ($qstart2, $qend2) = split (/=/, $range);
                my $match_flag = 0;
                if (($pre_qend2 == 0) and ($qstart2 > $min_align_len)){
                    $match_flag = 1;
                }
                elsif (($pre_qend2 > 0) and ($qstart2 - $pre_qend2 > $min_align_len)){
                    $match_flag = 2;
                }
                if ($match_flag > 0){
                    my %overlap;
                    my %overlap_rank;
                    foreach my $tname (keys %{$broken_query_seq_2{$qname}}){
                        next if ($tname eq $tname2);
                        my $overlap = 0;
                        foreach my $strand (keys %{${$broken_query_seq_2{$qname}}{$tname}}){
                            foreach my $qstart (sort {$a <=> $b} keys %{${${$broken_query_seq_2{$qname}}{$tname}}{$strand}}){
                                my $qend = ${${${$broken_query_seq_2{$qname}}{$tname}}{$strand}}{$qstart};
                                if ($match_flag == 1){
                                    last if ($qstart >= $qstart2);
                                    if (($qend <= $qstart2) and ($qend - $qstart + 1 >= $min_align_len)){
                                        $overlap += $qend - $qstart + 1;
                                        ${${$overlap{$tname}}{$strand}}{$qstart} = $qend;
                                    }
                                    elsif (($qend > $qstart2) and (($qstart2 - $qstart + 1 >= $min_align_len) or ($qstart2 - $qstart >= ($qend - $qstart) * 0.5))){
                                        $overlap += $qstart2 - $qstart + 1;
                                        ${${$overlap{$tname}}{$strand}}{$qstart} = $qend;
                                    }
                                }
                                elsif ($match_flag == 2){
                                    next if ($qend <= $pre_qend2);
                                    last if ($qstart >= $qstart2);
                                    if (($qstart >= $pre_qend2) and ($qend <= $qstart2) and ($qend - $qstart + 1 >= $min_align_len)){
                                        $overlap += $qend - $qstart + 1;
                                        ${${$overlap{$tname}}{$strand}}{$qstart} = $qend;
                                    }
                                    elsif (($pre_qend2 >= $qstart) and ($qend <= $qstart2) and ($qend - $pre_qend2 + 1 >= $min_align_len) and ($qend - $pre_qend2 >= ($qend - $qstart) * 0.5)){
                                        $overlap += $qend - $pre_qend2 + 1;
                                        ${${$overlap{$tname}}{$strand}}{$qstart} = $qend;
                                    }
                                    elsif (($pre_qend2 >= $qstart) and ($qend >= $qstart2) and ($qstart2 - $pre_qend2 + 1 >= $min_align_len) and ($qstart2 - $pre_qend2 >= ($qend - $qstart) * 0.5)){
                                        $overlap += $qstart2 - $pre_qend2 + 1;
                                        ${${$overlap{$tname}}{$strand}}{$qstart} = $qend;
                                    }
                                    elsif (($qstart >= $pre_qend2) and ($qend >= $qstart2) and (($qstart2 - $qstart + 1 >= $min_align_len) and ($qstart2 - $qstart >= ($qend - $qstart) * 0.5))){
                                        $overlap += $qstart2 - $qstart + 1;
                                        ${${$overlap{$tname}}{$strand}}{$qstart} = $qend;
                                    }
                                }
                            }
                            $overlap_rank{$overlap} = "$tname==$strand" if ($overlap > 0);
                        }
                    }
                    foreach my $value (sort {$b <=> $a} keys %overlap_rank){
                        my ($tname3, $strand3) = split (/==/, $overlap_rank{$value});
                        foreach my $qstart3 (sort {$a <=> $b} keys %{${$overlap{$tname3}}{$strand3}}){
                            my $qend3 = ${${$overlap{$tname3}}{$strand3}}{$qstart3};
                            if (!exists ${$trans_pos{$qname}}{$qstart2}){
                                ${$trans_pos{$qname}}{$qstart2} = "$tname3==$qstart3==$qend3==$strand3";
                            }
                            else{
                                ${$trans_pos{$qname}}{$qstart2} .= "||$tname3==$qstart3==$qend3==$strand3";
                            }
                        }
                        last;
                    }
                }
                $pre_qend2 = $qend2;
                $pre_qstart2 = $qstart2;
            }
            if ($qlen >= $pre_qend2 + $min_align_len){
                my %overlap;
                my %overlap_rank;
                foreach my $tname (keys %{$broken_query_seq_2{$qname}}){
                    next if ($tname eq $tname2);
                    my $overlap = 0;
                    foreach my $strand (keys %{${$broken_query_seq_2{$qname}}{$tname}}){
                        foreach my $qstart (sort {$a <=> $b} keys %{${${$broken_query_seq_2{$qname}}{$tname}}{$strand}}){
                            my $qend = ${${${$broken_query_seq_2{$qname}}{$tname}}{$strand}}{$qstart};
                            next if ($qend <= $pre_qend2);
                            if (($pre_qend2 <= $qstart) and ($qend - $qstart + 1 >= $min_align_len)){
                                $overlap += $qend - $qstart + 1;
                                ${${$overlap{$tname}}{$strand}}{$qstart} = $qend;
                            }
                            elsif (($pre_qend2 >= $qstart) and (($qend - $pre_qend2 + 1 >= $min_align_len) and ($qend - $pre_qend2 >= ($qend - $qstart) * 0.5))){
                                $overlap += $qend - $pre_qend2 + 1;
                                ${${$overlap{$tname}}{$strand}}{$qstart} = $qend;
                            }
                        }
                        $overlap_rank{$overlap} = "$tname==$strand" if ($overlap > 0);
                    }
                }
                foreach my $value (sort {$b <=> $a} keys %overlap_rank){
                    my ($tname3, $strand3) = split (/==/, $overlap_rank{$value});
                    foreach my $qstart3 (sort {$a <=> $b} keys %{${$overlap{$tname3}}{$strand3}}){
                        my $qend3 = ${${$overlap{$tname3}}{$strand3}}{$qstart3};
                        if (!exists ${$trans_pos{$qname}}{$pre_qend2}){
                            ${$trans_pos{$qname}}{$pre_qend2} = "$tname3==$qstart3==$qend3==$strand3";
                        }
                        else{
                            ${$trans_pos{$qname}}{$pre_qend2} .= "||$tname3==$qstart3==$qend3==$strand3";
                        }
                    }
                    last;
                }
            }
        }
        if (scalar keys %{${$broken_query_seq_2{$qname}}{$tname2}} > 1){
            my $strand3 = 'Plus';
            $strand3 = 'Minus' if ($strand2 eq 'Plus');
            my $pre_qend2 = 0;
            my $pre_qstart2 = 0;
            foreach my $range (@start_end){
                my ($qstart2, $qend2) = split (/=/, $range);
                my $match_flag = 0;
                if (($pre_qend2 == 0) and ($qstart2 > $min_align_len) and (!exists ${$trans_pos{$qname}}{$qstart2})){
                    $match_flag = 1;
                }
                elsif (($pre_qend2 > 0) and ($qstart2 - $pre_qend2 > $min_align_len) and (!exists ${$trans_pos{$qname}}{$qstart2})){
                    $match_flag = 2;
                }
                if ($match_flag > 0){
                    foreach my $qstart (sort {$a <=> $b} keys %{${${$broken_query_seq_2{$qname}}{$tname2}}{$strand3}}){
                        my $qend = ${${${$broken_query_seq_2{$qname}}{$tname2}}{$strand3}}{$qstart};
                        my $hit_flag = 0;
                        if ($match_flag == 1){
                            last if ($qstart >= $qstart2);
                            if (($qend <= $qstart2) and ($qend - $qstart + 1 >= $min_align_len)){
                                $hit_flag = 1;
                            }
                            elsif (($qend > $qstart2) and (($qstart2 - $qstart + 1 >= $min_align_len) or ($qstart2 - $qstart >= ($qend - $qstart) * 0.5))){
                                $hit_flag = 1;
                            }
                        }
                        elsif ($match_flag == 2){
                            next if ($qend <= $pre_qend2);
                            last if ($qstart >= $qstart2);
                            if (($qstart >= $pre_qend2) and ($qend <= $qstart2) and ($qend - $qstart + 1 >= $min_align_len)){
                                $hit_flag = 1;
                            }
                            elsif (($pre_qend2 >= $qstart) and ($qend <= $qstart2) and ($qend - $pre_qend2 + 1 >= $min_align_len) and ($qend - $pre_qend2 >= ($qend - $qstart) * 0.5)){
                                $hit_flag = 1;
                            }
                            elsif (($pre_qend2 >= $qstart) and ($qend >= $qstart2) and ($qstart2 - $pre_qend2 + 1 >= $min_align_len) and ($qstart2 - $pre_qend2 >= ($qend - $qstart) * 0.5)){
                                $hit_flag = 1;
                            }
                            elsif (($qstart >= $pre_qend2) and ($qend >= $qstart2) and (($qstart2 - $qstart + 1 >= $min_align_len) and ($qstart2 - $qstart >= ($qend - $qstart) * 0.5))){
                                $hit_flag = 1;
                            }
                        }
                        if ($hit_flag == 1){
                            if (!exists ${$inv_pos{$qname}}{$qstart2}){
                                ${$inv_pos{$qname}}{$qstart2} = "$strand3==$qstart==$qend";
                            }
                            else{
                                ${$inv_pos{$qname}}{$qstart2} .= "||$strand3==$qstart==$qend";
                            }
                        }
                    }
                }
                $pre_qend2 = $qend2;
                $pre_qstart2 = $qstart2;
            }
            if (($qlen >= $pre_qend2 + $min_align_len) and (!exists ${$trans_pos{$qname}}{$pre_qend2})){
                my %overlap;
                my %overlap_rank;
                foreach my $qstart (sort {$a <=> $b} keys %{${${$broken_query_seq_2{$qname}}{$tname2}}{$strand3}}){
                    my $qend = ${${${$broken_query_seq_2{$qname}}{$tname2}}{$strand3}}{$qstart};
                    my $hit_flag = 0;
                    next if ($qend <= $pre_qend2);
                    if (($pre_qend2 <= $qstart) and ($qend - $qstart + 1 >= $min_align_len)){
                        $hit_flag = 1;
                    }
                    elsif (($pre_qend2 >= $qstart) and (($qend - $pre_qend2 + 1 >= $min_align_len) and ($qend - $pre_qend2 >= ($qend - $qstart) * 0.5))){
                        $hit_flag = 1;
                    }
                    if ($hit_flag == 1){
                        if (!exists ${$inv_pos{$qname}}{$pre_qend2}){
                            ${$inv_pos{$qname}}{$pre_qend2} = "$strand3==$qstart==$qend";
                        }
                        else{
                            ${$inv_pos{$qname}}{$pre_qend2} .= "||$strand3==$qstart==$qend";
                        }
                    }
                }
            }
        }
    }
}


my $out_stat_file = $prefix_out . '.stat.txt';
my $corrected_contig_file = $prefix_out . '.corrected.fa';
my $out_misassemble_query_list = $prefix_out . '.misassemble.list.txt';

my $out_truelyassemble_query_list = $prefix_out . '.true-assemble.list.txt';

my $total_misassemble_num = 0;
my $misassemble_query_num = 0;
my @misassemble_query_name;
my $total_query_num = scalar keys %query_seq;
my $truely_assemble_query_num = scalar keys %assemble_query;

my $corrected_total_len = 0;
my @corrected_seq_len;
my @misassembly_list;
open (OUT, "> $corrected_contig_file") if ($error_correct == 1);
open (OUT2, "> $out_misassemble_query_list");
open (OUT3, "> $out_truelyassemble_query_list");
foreach my $qname (keys %query_seq){
    if (exists $assemble_query{$qname} and !exists $broken_query_seq_2{$qname}){
        if ($query_length{$qname} >= $min_align_len){
            if ($error_correct == 1){
                my $qseq = $query_seq{$qname};
                my ($qstart, $qend) = split (/=/, $assemble_query{$qname});
                if ($qend < $query_length{$qname} - 2) {
                    substr ($qseq, $qend, $query_length{$qname} - $qend, '');
                }
                if ($qstart > 3){
                    substr ($qseq, 0, $qstart - 1, '');
                }
                print OUT ">$qname\n";
                print OUT $qseq, "\n";
            }
            print OUT3 "$qname\n";
        }
    }
    if (!exists $assemble_query{$qname} and !exists $broken_query_seq_2{$qname} and exists $misassemble_query{$qname}){
        $misassemble_query_num ++;
        $total_misassemble_num ++;
        $small_cov ++;
        push @misassembly_list, $qname;
        my $pre_qstart = 0;
        my $qlen = length $query_seq{$qname};
        foreach my $qstart (sort {$a <=> $b} keys %{$misassemble_query{$qname}}){
            $qstart = 1 if ($qstart < 1);
            my $broken_len = ${$misassemble_query{$qname}}{$qstart};
            if ($pre_qstart == 0){
                if ($qstart > 1){
                    if ($qstart - 1 >= $min_align_len){
                        $corrected_total_len += $qstart - 1;
                        push @corrected_seq_len, $qstart - 1;
                    }
                    $pre_qstart = $qstart;
                    next;
                }
            }
            if ($qstart - $pre_qstart >= $min_align_len){
                $corrected_total_len += $qstart - $pre_qstart;
                push @corrected_seq_len, $qstart - $pre_qstart;
            }
            $pre_qstart = $qstart;
        }
        if (($pre_qstart < $qlen) and ($qlen - $pre_qstart + 1 >= $min_align_len)){
            $corrected_total_len += $qlen - $pre_qstart + 1;
            push @corrected_seq_len, $qlen - $pre_qstart + 1;
        }
    }
}
close (OUT3);
foreach my $qname (keys %broken_query_seq_2){
    next if (exists $assemble_query{$qname});
    my $qlen = length $query_seq{$qname};
    $misassemble_query_num ++;
    if (exists $broken_query_cov{$qname}){
        my ($tname2, $strand2, $overlap2) = split (/==/, $broken_query_cov{$qname});
        my $pre_qend = 0;
        my $pre_tend = 0;
        my $pre_tstart = 0;
        my $sub_num = 0;
        my $broken_num = 0;
        foreach my $qstart (sort {$a <=> $b} keys %{${${$broken_query_seq_2{$qname}}{$tname2}}{$strand2}}){
            $qstart = 1 if ($qstart < 1);
            my $qend = ${${${$broken_query_seq_2{$qname}}{$tname2}}{$strand2}}{$qstart};
            $qend = $qlen if ($qend > $qlen);
            next if ($pre_qend >= $qend);
            $broken_num ++;
            $sub_num ++;
            my ($tstart, $tend) = (0, 0);
            ($tstart, $tend) = split (/=/, ${${${$broken_target_seq{$qname}}{$tname2}}{$strand2}}{$qstart}) if (exists ${${${$broken_target_seq{$qname}}{$tname2}}{$strand2}}{$qstart});
            if ($pre_qend > 0){
                if (($pre_tstart > 0) and ($tstart > 0)){
                    if ($strand2 eq 'plus'){
                        if ((abs ($pre_tend - $tstart) >= 1000) or (abs ($pre_qend - $qstart) >= 1000)){
                            $local_error_L ++;
                        }
                    }
                    elsif ($strand2 eq 'Minus'){
                        if ((abs ($pre_tstart - $tend) >= 1000) or (abs ($pre_qend - $qstart) >= 1000)){
                            $local_error_L ++;
                        }
                    }
                }
                else{
                    if (abs ($pre_qend - $qstart) >= 1000){
                        $local_error_L ++;
                    }
                }
            }
            if ($qend - $qstart + 1 >= $min_align_len){
                print OUT2 "$qname\t$qstart\t$qend\t$query_length{$qname}\t$tname2\tlocal-error\n";
                $corrected_total_len += $qend - $qstart + 1;
                push @corrected_seq_len, $qend - $qstart + 1;
            }
            if ($error_correct == 1){
                if ($qstart <= $pre_qend){
                    $qstart = $pre_qend + 1;
                }
                if ($qend - $qstart + 1 >= $min_align_len){
                    print OUT ">$qname-$sub_num\n";
                    print OUT substr ($query_seq{$qname}, $qstart - 1, $qend - $qstart + 1), "\n";
                }
            }
            $pre_qend = $qend;
            $pre_tend = $tend;
            $pre_tstart = $tstart;
        }
        $total_misassemble_num ++ if ($broken_num == 1);
        $total_misassemble_num += $broken_num - 1 if ($broken_num >= 2);
        $small_cov ++ if ($broken_num == 1);
        $local_error += $broken_num - 1 if ($broken_num >= 2); 
    }
    elsif (exists $broken_query_uncov{$qname}){
        my ($tname2, $strand2, $overlap2) = split (/==/, $broken_query_uncov{$qname});
        my $pre_qend = 0;
        my $pre_tend = 0;
        my $pre_tstart = 0;
        my $sub_num = 0;
        my $broken_num = 0;
        my $troc_num = 0;
        my $inv_num = 0;
        foreach my $qstart (sort {$a <=> $b} keys %{${${$broken_query_seq_2{$qname}}{$tname2}}{$strand2}}){
            $qstart = 1 if ($qstart < 1);
            my $qend = ${${${$broken_query_seq_2{$qname}}{$tname2}}{$strand2}}{$qstart};
            $qend = $qlen if ($qend > $qlen);
            next if ($pre_qend >= $qend);
            my ($tstart, $tend) = (0, 0);
            $broken_num ++;
            if ((exists ${$trans_pos{$qname}}{$qstart}) or (exists ${$trans_pos{$qname}}{$qend})){
                $translocate ++;
                $troc_num ++;
                my @align = ();
                if (exists ${$trans_pos{$qname}}{$qstart}){
                    @align = split (/\|\|/, ${$trans_pos{$qname}}{$qstart});
                }
                else{
                    @align = split (/\|\|/, ${$trans_pos{$qname}}{$qend});
                }
                my $pre_qend3 = 0;
                my $pre_tend3 = 0;
                my $pre_tstart3 = 0;
                foreach my $align (@align){
                    $broken_num ++;
                    my ($tname3, $qstart3, $qend3, $strand3) = split (/==/, $align);
                    my ($tstart3, $tend3) = (0, 0);
                    ($tstart3, $tend3) = split (/=/, ${${${$broken_target_seq{$qname}}{$tname3}}{$strand3}}{$qstart3}) if (exists ${${${$broken_target_seq{$qname}}{$tname3}}{$strand3}}{$qstart3});
                    if ($pre_qend3 > 0){
                        if (($pre_tstart3 > 0) and ($tstart3 > 0)){
                            if ($strand2 eq 'plus'){
                                if ((abs ($pre_tend3 - $tstart3) >= 1000) or (abs ($pre_qend3 - $qstart3) >= 1000)){
                                    $local_error_L ++;
                                }
                            }
                            elsif ($strand2 eq 'Minus'){
                                if ((abs ($pre_tstart3 - $tend3) >= 1000) or (abs ($pre_qend3 - $qstart3) >= 1000)){
                                    $local_error_L ++;
                                }
                            }
                        }
                        else{
                            if (abs ($pre_qend3 - $qstart3) >= 1000){
                                $local_error_L ++;
                            }
                        }
                    }
                    if ($pre_qend3 == 0){
                        if ($qend3 - $qstart3 + 1 >= $min_align_len){
                            print OUT2 "$qname\t$qstart3\t$qend3\t$query_length{$qname}\t$tname3\ttranslocation\n";
                            $corrected_total_len += $qend - $qstart + 1;
                            push @corrected_seq_len, $qend - $qstart + 1;
                        }
                    }
                    else{
                        if ($qend3 - $qstart3 + 1 >= $min_align_len){
                            print OUT2 "$qname\t$qstart3\t$qend3\t$query_length{$qname}\t$tname3\tlocal-misassembly\n";
                            $corrected_total_len += $qend - $qstart + 1;
                            push @corrected_seq_len, $qend - $qstart + 1;
                        }
                    }
                    $sub_num ++;
                    if ($error_correct == 1){
                        if ($qend3 - $qstart3 + 1 >= $min_align_len){
                            print OUT ">$qname-$sub_num\n";
                            print OUT substr ($query_seq{$qname}, $qstart3 - 1, $qend3 - $qstart3 + 1), "\n";
                        }
                    }
                    $pre_qend3 = $qend3;
                    $pre_tstart3 = $tstart3;
                    $pre_tend3 = $tend3;
                }
            }
            elsif ((exists ${$inv_pos{$qname}}{$qstart}) or (exists ${$inv_pos{$qname}}{$qend})){
                $inversion ++;
                $inv_num ++;
                my @align = ();
                if (exists ${$inv_pos{$qname}}{$qstart}){
                    @align = split (/\|\|/, ${$inv_pos{$qname}}{$qstart});
                }
                else{
                    @align = split (/\|\|/, ${$inv_pos{$qname}}{$qend});
                }
                my $pre_qend3 = 0;
                my $pre_tend3 = 0;
                my $pre_tstart3 = 0;
                foreach my $align (@align){
                    $broken_num ++;
                    my ($strand3, $qstart3, $qend3) = split (/==/, $align);
                    my ($tstart3, $tend3) = (0, 0);
                    ($tstart3, $tend3) = split (/=/, ${${${$broken_target_seq{$qname}}{$tname2}}{$strand3}}{$qstart3}) if (exists ${${${$broken_target_seq{$qname}}{$tname2}}{$strand3}}{$qstart3});
                    if ($pre_qend3 > 0){
                        if (($pre_tstart3 > 0) and ($tstart3 > 0)){
                            if ($strand2 eq 'plus'){
                                if ((abs ($pre_tend3 - $tstart3) >= 1000) or (abs ($pre_qend3 - $qstart3) >= 1000)){
                                    $local_error_L ++;
                                }
                            }
                            elsif ($strand2 eq 'Minus'){
                                if ((abs ($pre_tstart3 - $tend3) >= 1000) or (abs ($pre_qend3 - $qstart3) >= 1000)){
                                    $local_error_L ++;
                                }
                            }
                        }
                        else{
                            if (abs ($pre_qend3 - $qstart3) >= 1000){
                                $local_error_L ++;
                            }
                        }
                    }
                    if ($pre_qend3 == 0){
                        if ($qend - $qstart + 1 >= $min_align_len){
                            print OUT2 "$qname\t$qstart3\t$qend3\t$query_length{$qname}\t$tname2\tinversion\n";
                            $corrected_total_len += $qend - $qstart + 1;
                            push @corrected_seq_len, $qend - $qstart + 1;
                        }
                    }
                    else{
                        if ($qend - $qstart + 1 >= $min_align_len){
                            print OUT2 "$qname\t$qstart3\t$qend3\t$query_length{$qname}\t$tname2\tlocal-misassembly\n";
                            $corrected_total_len += $qend - $qstart + 1;
                            push @corrected_seq_len, $qend - $qstart + 1;
                        }
                    }
                    $sub_num ++;
                    if ($error_correct == 1){
                        if ($qend3 - $qstart3 + 1 >= $min_align_len){
                            print OUT ">$qname-$sub_num\n";
                            print OUT substr ($query_seq{$qname}, $qstart3 - 1, $qend3 - $qstart3 + 1), "\n";
                        }
                    }
                    $pre_qend3 = $qend3;
                    $pre_tstart3 = $tstart3;
                    $pre_tend3 = $tend3;
                }
            }
            else{
                ($tstart, $tend) = split (/=/, ${${${$broken_target_seq{$qname}}{$tname2}}{$strand2}}{$qstart}) if (exists ${${${$broken_target_seq{$qname}}{$tname2}}{$strand2}}{$qstart});
                if ($pre_qend > 0){
                    if (($pre_tstart > 0) and ($tstart > 0)){
                        if ($strand2 eq 'plus'){
                            if ((abs ($pre_tend - $tstart) >= 1000) or (abs ($pre_qend - $qstart) >= 1000)){
                                $local_error_L ++;
                            }
                        }
                        elsif ($strand2 eq 'Minus'){
                            if ((abs ($pre_tstart - $tend) >= 1000) or (abs ($pre_qend - $qstart) >= 1000)){
                                $local_error_L ++;
                            }
                        }
                    }
                    else{
                        if (abs ($pre_qend - $qstart) >= 1000){
                            $local_error_L ++;
                        }
                    }
                }
            }
            if ($qend - $qstart + 1 >= $min_align_len){
                print OUT2 "$qname\t$qstart\t$qend\t$query_length{$qname}\t$tname2\tlocal-misassembly\n";
                $corrected_total_len += $qend - $qstart + 1;
                push @corrected_seq_len, $qend - $qstart + 1;
            }
            $sub_num ++;
            if ($error_correct == 1){
                if ($qstart <= $pre_qend){
                    $qstart = $pre_qend + 1;
                }
                if ($qend - $qstart + 1 >= $min_align_len){
                    print OUT ">$qname-$sub_num\n";
                    print OUT substr ($query_seq{$qname}, $qstart - 1, $qend - $qstart + 1), "\n";
                }
            }
            $pre_qend = $qend;
            $pre_tend = $tend;
            $pre_tstart = $tstart;
        }
        $total_misassemble_num ++ if ($broken_num == 1);
        $total_misassemble_num += $broken_num - 1 if ($broken_num >= 2);
        $small_cov ++ if ($broken_num == 1);
        my $local_num = $broken_num - $troc_num - $inv_num;
        $local_error += $local_num - 1;
    }
}
close (OUT) if ($error_correct == 1);
print OUT2 "\n<<Query contigs of which the total length of the alignment overlaps is smaller than min_contig_length>>\n" if (@misassembly_list > 0);
foreach (@misassembly_list){
    print OUT2 "$_\n";
}
close (OUT2);



my $total_len = 0;                      # calculate N50 and total assembly length of the assemblies
my @seq_len;
my $count_len = 0;
my $half_len = 0;
my $N50 = 0;

foreach my $qname (keys %query_length){
    if ($query_length{$qname} >= $min_align_len){
        $total_len += $query_length{$qname};
        push @seq_len, $query_length{$qname};
        if (exists $assemble_query{$qname} and !exists $broken_query_seq_2{$qname}){
            $corrected_total_len += $query_length{$qname};
            push @corrected_seq_len, $query_length{$qname};
        }
        elsif (!exists $assemble_query{$qname} and !exists $broken_query_seq_2{$qname} and !exists $misassemble_query{$qname}){
            $corrected_total_len += $query_length{$qname};
            push @corrected_seq_len, $query_length{$qname};
        }
    }
}
@seq_len = sort {$b<=>$a} @seq_len;
for (my $i = 0; $i < @seq_len; $i++){
    $count_len += $seq_len[$i];
    if (($count_len >= $total_len / 2) and ($half_len == 0)){
        $N50 = $seq_len[$i];
        $half_len = $seq_len[$i];
    }
}

$count_len = 0;                         # calculate N50 and total assembly length of the assemblies that are broken at mis-assemble sites
$half_len = 0;
my $corrected_N50 = 0;
@corrected_seq_len = sort {$b<=>$a} @corrected_seq_len;
for (my $i = 0; $i < @corrected_seq_len; $i++){
    $count_len += $corrected_seq_len[$i];
    if (($count_len >= $corrected_total_len / 2) and ($half_len == 0)){
        $corrected_N50 = $corrected_seq_len[$i];
        $half_len = $corrected_seq_len[$i];
    }
}

open (OUT, "> $out_stat_file");
    my $unmapped_query_num = $total_query_num - $truely_assemble_query_num - $misassemble_query_num;
    my $rate_misassemble_query = 0;
    my $rate_truely_assemble_query = 0;
    $rate_misassemble_query = int ($misassemble_query_num / ($truely_assemble_query_num + $misassemble_query_num) * 1000) / 10 if ($truely_assemble_query_num + $misassemble_query_num > 0);
    $rate_truely_assemble_query = int ($truely_assemble_query_num / ($truely_assemble_query_num + $misassemble_query_num) * 1000) / 10 if ($truely_assemble_query_num + $misassemble_query_num > 0);
    my $rate_unmapped_query = int ($unmapped_query_num / $total_query_num * 1000) / 10;
    my $local_error_S = $local_error - $local_error_L;
    print OUT "Number of contigs: $total_query_num\n";
    print OUT "Total number of misassembled contigs = $misassemble_query_num ($rate_misassemble_query%)\n";
    print OUT "Total number of truely assembled contigs = $truely_assemble_query_num ($rate_truely_assemble_query%)\n";
    print OUT "Total number of misassembled events = $total_misassemble_num\n";
    print OUT "Number of unmapped contigs = $unmapped_query_num ($rate_unmapped_query%)\n\n";
    
    print OUT "Misassembly events on the same chromosome (break point distance > $max_indel_size bp, < 1000 bp): $local_error_S\n";
    print OUT "Misassembly events (Relocation) on the same chromosome (break point distance >= 1000 bp): $local_error_L\n";
    print OUT "Misassembly events (Translocation): $translocate\n";
    print OUT "Misassembly events (Inversion): $inversion\n";
    print OUT "Misassembly events (alignment coverage < $min_coverage%): $small_cov\n";
    print OUT "Medium indels (> 5 bp and <= $max_indel_size bp): $middle_indel\n\n";
    
    print OUT "Total contig length: $total_len\n";
    print OUT "Contig N50: $N50\n";
    print OUT "Corrected contig N50: $corrected_N50\n";
close (OUT);

print STDERR "Total number of contigs = $total_query_num\n";
print STDERR "Total number of misassembled contigs = $misassemble_query_num ($rate_misassemble_query%)\n";
print STDERR "Total number of misassembled events = $total_misassemble_num\n";
print STDERR "Number of unmapped contigs = $unmapped_query_num ($rate_unmapped_query%)\n";
print STDERR "N50: $N50\tCorrected N50: $corrected_N50\n";
