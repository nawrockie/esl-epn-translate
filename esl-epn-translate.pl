#!/usr/bin/env perl
# 
# esl-test-cds-translate-vs-fetch.pl: 
#
#  Given a sequence file output from esl-fetch-cds.pl fetch the AA
#  sequence for each CDS and test that the CDS sequence translates to
#  exactly that sequence.
#
#  It's important that the sequence file was created by
#  esl-fetch-cds.pl because it relies on some naming conventions of
#  that script.  In particular, the first ':' delimited token of the
#  sequence name should be the protein accession.
#
# This script uses BioEasel's SqFile module and the 'idfetch' program to 
# fetch the protein sequence (hard-coded path).
# 
# EPN, Tue Mar 31 12:50:15 2015

use strict;
use Getopt::Long;
use Bio::Easel::SqFile;

my $in_fafile        = "";  # name of input fasta file
my $idfetch          = "/netopt/ncbi_tools64/bin/idfetch";
my $nall             = 11; # number of total tests
my $nsubset          = 9;  # number of tests run if -subset is used
my $do_all           = 0;  # changed to '1' with -a to print all sequences, not just those that fail >= 1 test or have ambig chars
my $do_verbose       = 0;  # changed to '1' with -v, output fetched and translated protein sequences
my $skip_incompletes = 0;  # changed to '1' with -skipinc, output fetched and translated protein sequences
my $do_subset        = 0;  # changed to '1' with -subset, run only first $nsubset tests instead of all $nall
my $do_compare_input = 0;  # changed to '1' with -incompare, input sequences were created by dnaorg_compare_genomes.pl

&GetOptions( "v"         => \$do_verbose, 
             "a"         => \$do_all, 
             "incompare" => \$do_compare_input,
             "subset"    => \$do_subset,
             "skipinc"   => \$skip_incompletes) || die "ERROR unknown option";

my $usage;
$usage  = "esl-test-cds-translate-vs-fetch.pl [OPTIONS] <input fasta file output from esl-fetch-cds.pl>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-v         : be verbose; output translated and fetched protein sequences\n";
$usage .= "\t\t-a         : print all sequences, even those that pass all tests and have 0 ambig chars\n";
$usage .= "\t\t-subset    : run only the first $nsubset tests, not all $nall\n";
$usage .= "\t\t-incompare : input file was created by dnaorg_compare_genomes.pl -protid -codonstart\n";
$usage .= "\t\t-skipinc   : skip examination of incomplete CDS'\n";

if(scalar(@ARGV) != 1) { die $usage; }
($in_fafile) = @ARGV;

if(! -s $in_fafile) { die "ERROR, $in_fafile does not exist"; }

# open sequence file, index it, and determine the number of sequences in it
my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $in_fafile });
my $nseq = $sqfile->nseq_ssi;
my @testdesc_A = (); # list of descriptions for each test

# print headers for tabular output
my @complete_nfail_A  = (); # [0..$j..$last_test-1] number of failures for complete CDS for test $j+1
my @incomplete_nfail_A= (); # [0..$j..$last_test-1] number of failures for incomplete CDS for test $j+1
my $ncomplete        = 0; # number of complete CDS, total
my $nincomplete      = 0; # number of incomplete CDS, total
my $ncomplete_fail   = 0; # number of complete CDS that had at least 1 failure
my $nincomplete_fail = 0; # number of incomplete CDS that had at least 1 failure
my $last_test = ($do_subset) ? $nsubset : $nall;
my $paccn_w   = length("protein-accession");
my $ntaccn_w  = length("nt-accession");

#$paccn_w  = 22;
#$ntaccn_w = 22;

# for each sequence in $in_fafile:
# 1. fetch the CDS sequence from $in_fafile
# 2. determine if CDS sequence is incomplete on 5' and/or 3' end
# 3. translate the CDS sequence
# 4. fetch the protein sequence using idfetch
#
# Then test the translated CDS and fetched proteins:

my @c_toprint_A       = (); # array of output lines for complete CDS
my @icstart_toprint_A = (); # array of output lines for incomplete CDS, incomplete start, complete stop
my @icstop_toprint_A  = (); # array of output lines for incomplete CDS, complete start, incomplete stop
my @icboth_toprint_A  = (); # array of output lines for incomplete CDS, incomplete start, incomplete stop
my $toprint;

# initialize
for(my $i = 0; $i < $last_test; $i++) { 
  $complete_nfail_A[$i]   = 0;
  $incomplete_nfail_A[$i] = 0;
}

# write the test descriptions:
push(@testdesc_A, "Length of CDS (DNA) is a multiple of 3 (or stop is incomplete)"); 
push(@testdesc_A, "Fetched protein starts with an M (or is annot. as incomplete on 5' end)"); 
push(@testdesc_A, "CDS ends with a stop codon (or is annot. as incomplete on 3' end)"); 
push(@testdesc_A, "CDS has 0 Ns that make AA assignment ambiguous"); 
push(@testdesc_A, "Translated CDS and fetched protein are identical length\n#              (ignoring final codon, if it's a stop or incomplete)"); 
push(@testdesc_A, "Translated CDS and fetched protein are identical sequence\n#              (ignoring hangovers if lengths differ and final codon (if stop) and codons with ambiguous nts)"); 
push(@testdesc_A, "No internal stops in translated CDS"); 
push(@testdesc_A, "No internal stops in fetched protein"); 
push(@testdesc_A, "CDS has a linked protein accession"); 
push(@testdesc_A, "CDS is complete on the 5' end (not annotated as incomplete on 5' end)");
push(@testdesc_A, "CDS is complete on the 3' end (not annotated as incomplete on 3' end)");

# for each sequence in $in_fafile:
for(my $i = 0; $i < $nseq; $i++) { 
  # 1. fetch the CDS sequence from $in_fafile
  my ($cds_name) = $sqfile->fetch_seq_name_given_ssi_number($i);

  my ($cds_name2, $cds_seq) = split(/\n/, $sqfile->fetch_seq_to_fasta_string($cds_name, -1));
  $cds_name2 =~ s/^\>//;
  my $cds_desc = $cds_name2; 
  $cds_desc =~ s/^\S+//;
  $cds_desc =~ s/^\s+//;
  $cds_name2 =~ s/\s+.+$//;

  # default input example:
  # >AAC62262.1:codon_start1:AC005031.1:55514:55564:-:AC005031.1:61328:61438:-

  # -incompare input:
  # >NC_004004 NC_004004:1059:8027:+: product:pol protein protein_id:ref|NP_658990.1| codon_start:none-annotated
  # >NC_009942 NC_009942:97:3552:+:NC_009942:3552:3683:+: product:truncated polyprotein protein_id:ref|YP_006485882.1| codon_start:none-annotated

  # sanity check
  if($cds_name ne $cds_name2) { die "ERROR, unexpected error, name mismatch ($cds_name ne $cds_name2)"; }

  # 1. isolate the protein accession and codon_start value
  my $paccn;
  my $codon_start;
  if($do_compare_input) { 
    $paccn = $cds_desc;
    # NC_004004:1059:8027:+: product:pol protein protein_id:ref|NP_658990.1| codon_start:none-annotated
    $paccn =~ s/^.+protein\_id\://;
    # ref|NP_658990.1| codon_start:none-annotated
    $paccn =~ s/\s+.+$//;
    # ref|NP_658990.1|
    $paccn =~ s/^\w+\|//;
    # NP_658990.1|
    $paccn =~ s/\|.*$//;
    # NP_658990.1

    $codon_start = $cds_desc;
    # NC_004004:1059:8027:+: product:pol protein protein_id:ref|NP_658990.1| codon_start:none-annotated
    $codon_start =~ s/^.+codon\_start\://;
    # none-annotated
    $codon_start =~ s/\s+.+$//;
    # none-annotated
    if($codon_start eq "none-annotated" || $codon_start eq "-") { 
      $codon_start = 1;
    }
    elsif($codon_start ne "1" && $codon_start ne "2" && $codon_start ne "3") {
      die "ERROR (-incompare) unable to parse codon_start in desc: $cds_desc (is -incompare appropriate to use for this file?)";
    }
  }
  else { # default case 
    $paccn = $cds_name;
    # >AAC62262.1:codon_start1:AC005031.1:55514:55564:-:AC005031.1:61328:61438:-
    $paccn =~ s/^\>//;
    # AAC62262.1:codon_start1:AC005031.1:55514:55564:-:AC005031.1:61328:61438:-
    $paccn =~ s/\:.+$//;
    # AAC62262.1
    
    $codon_start = $cds_name;
    # >AAC62262.1:codon_start1:AC005031.1:55514:55564:-:AC005031.1:61328:61438:-
    $codon_start =~ s/^.+\:codon\_start//;
    # 1:AC005031.1:55514:55564:-:AC005031.1:61328:61438:-
    $codon_start =~ s/\:.+$//;
    # 1
    if($codon_start ne "1" && $codon_start ne "2" && $codon_start ne "3") {
      die "ERROR (default) unable to parse codon_start in $cds_name (codon_start: $codon_start, should you be using -incompare for this file?))";
    }
  }

  # 2. determine if CDS sequence is incomplete on 5' and/or 3' end
  my $cds_name_or_desc = ($do_compare_input) ? $cds_desc : $cds_name;
  my ($accn_coord_strand_str, $ic_start, $ic_stop) = characterizeCDS($cds_name_or_desc, $do_compare_input);
  my ($ntaccn, $mincoord, $maxcoord, $nexon, $strand2print) = split(":::", $accn_coord_strand_str);

  # 3. translate the CDS sequence
  my ($prot_translated, $n_N, $n_nonACGTN) = translateDNA($cds_seq, $codon_start);

  # 4. fetch the protein sequence using idfetch, if we have an accession that is (rarely, we won't and $paccn eq "-")
  my @fail_A  = (); # '0' if test x failed, else 0
  my $cds_len = length($cds_seq);
  my $cds_len_to_translate = $cds_len - ($codon_start - 1);
  my $errmsg = "";
  my $prot_fetched = "";

  if($paccn ne "-") { 
    # remove 'version' from $paccn
    my $prot_acconly = $paccn;
    $prot_acconly =~ s/\.\d+$//;
    # need to create a temporary file with the accession for idfetch
    my $tmp_acc_file = "tmp.$prot_acconly.acc";
    open(OUT, ">" . $tmp_acc_file);
    print OUT $prot_acconly . "\n";
    # fetch the sequence
    my $idfetch_output = `$idfetch -t 5 -c 1 -G $tmp_acc_file`;
    # remove name and newlines
    my @idfetch_output_A = split(/\n/, $idfetch_output);
    my $nlines = scalar(@idfetch_output_A);
    $prot_fetched = "";
    for(my $i = 1; $i < $nlines; $i++) { # skip first
      if($idfetch_output_A[$i] !~ m/^\>/ && 
         $idfetch_output_A[$i] =~ m/\w/) { 
        $prot_fetched .= $idfetch_output_A[$i]
      }
    }
    # remove temp file with accession
    unlink $tmp_acc_file;
  }    

  # Run the tests:
  # Test the translated CDS and fetched protein
  my $translated_len = length($prot_translated);
  my $fetched_len    = length($prot_fetched);
  my @translated_A   = split("", $prot_translated); 
  my @fetched_A      = split("", $prot_fetched);
    
  # Test: Is CDS length divisible by 3? 
  if($ic_stop) { 
    push(@fail_A, 0); # pass
  }
  else { 
    if(($cds_len_to_translate % 3) != 0) { push(@fail_A, 1); } # fail
    else                                 { push(@fail_A, 0); } # pass
  }
  
  # Test: Does fetched protein start with an M (if we have a protein and it's supposed to be complete at start)?
  if($ic_start || $paccn eq "-") { 
    push(@fail_A, 0); # pass 
  }
  else { 
    if($fetched_A[0] ne "M") { push(@fail_A, 1); } # fail
    else                     { push(@fail_A, 0); } # pass
  }
  
  # Test: Does translated CDS end with a stop codon?
  if($ic_stop) { 
    push(@fail_A, 0); # pass 
  }
  else { 
    if($translated_A[($translated_len-1)] ne "*") { push(@fail_A, 1); } # fail
    else                                          { push(@fail_A, 0); } # pass
  }
  
  # Test: Are there any nucleotide ambiguities that introduce an AA ambiguity?
  my $nambig_aa = 0;
  for(my $p = 0; $p < $translated_len; $p++) { 
    if($translated_A[$p] eq "?") { 
      $nambig_aa++;
    }
  }
  if($nambig_aa > 0) { push(@fail_A, 1); } # fail
  else               { push(@fail_A, 0); } # pass
  
  # Test: Do the protein lengths match?
  if($paccn ne "-") { 
    if(($translated_A[($translated_len-1)] eq "*") || # CDS ended with a stop codon
       ($translated_A[($translated_len-1)] eq "~")) { # CDS ended with a partial codon
      if(($translated_len-1) != $fetched_len) { push(@fail_A, 1); } # fail
      else                                    { push(@fail_A, 0); } # pass
    }
    else { # CDS does not end with a stop codon or partial codon
      if($translated_len != $fetched_len) { push(@fail_A, 1); } # fail
      else                                { push(@fail_A, 0); } # pass
    }
  } 
  else { # we don't have a fetched protein, so we 'pass' this test
    push(@fail_A, 0); # pass
  }  

  # Test: Are there any mismatches between translated CDS and fetched protein?
  my $min_len      = ($fetched_len < $translated_len) ? $fetched_len : $translated_len;
  my $nmismatch       = 0;
  my $nambig_mismatch = 0;
  for(my $p = 0; $p < $min_len; $p++) { 
    if($translated_A[$p] ne $fetched_A[$p]) { 
      if(($p == 0) && ($fetched_A[$p] eq "M")) { # we allow this
        ;
      }
      elsif(($translated_A[$p] eq "*") && ($fetched_A[$p] eq "X")) { # we allow this
        ;
      }
      elsif($translated_A[$p] eq "?") { 
        $nambig_mismatch++;
      }
      else { 
        $nmismatch++;
        if($errmsg ne "") { $errmsg .= "; "; }
        $errmsg .= sprintf("position %d mismatch %s ne %s (translated ne fetched)", $p+1, $translated_A[$p], $fetched_A[$p]);
        if($do_verbose) { print $paccn . " " . $errmsg . "\n"; }
      }
    }
  }
  if($paccn ne "-" && $nmismatch > 0) { push(@fail_A, 1); } # fail
  else                                { push(@fail_A, 0); } # pass

  # Test: Are there any internal stop codons in the translated CDS?
  my $nintstop_translated = 0;
  for(my $p = 0; $p < ($translated_len-1); $p++) { 
    if($translated_A[$p] eq "*") { 
      $nintstop_translated++;
    }
  }
  if($nintstop_translated > 0) { push(@fail_A, 1); } # fail
  else                         { push(@fail_A, 0); } # pass

  # Test: Are there any internal stop codons in the fetched protein?
  my $nintstop_fetched = 0;
  for(my $p = 0; $p < ($fetched_len-1); $p++) { 
    if($fetched_A[$p] eq "*") { 
      $nintstop_fetched++;
    }
  }
  if($paccn ne "-" && $nintstop_fetched > 0) { push(@fail_A, 1); } # fail
  else                                       { push(@fail_A, 0); } # pass

  # Test: do we have a valid protein accession? 
  if($paccn eq "-") { push(@fail_A, 1); } # fail
  else              { push(@fail_A, 0); } # pass

  # Test: are we complete on 5' end?
  if($ic_start) { push(@fail_A, 1); } # fail
  else          { push(@fail_A, 0); } # pass

  # Test: are we complete on 3' end?
  if($ic_stop) { push(@fail_A, 1); } # fail
  else         { push(@fail_A, 0); } # pass

  # store output for to print later
  my $any_fails = 0;
  for(my $j = 0; $j < $last_test; $j++) { 
    if($fail_A[$j]) { 
      $any_fails = 1; 
    }
  }
  $toprint = "";
  my $incomplete = "no";
  if($ic_start && $ic_stop) { $incomplete = "yes(both)"; }
  elsif($ic_start)          { $incomplete = "yes(start)"; }
  elsif($ic_stop)           { $incomplete = "yes(stop)"; }
  if($incomplete eq "no") { $ncomplete++;   }
  else                    { $nincomplete++; }

  if($any_fails) { 
    if($incomplete eq "no") { $ncomplete_fail++;   }
    else                    { $nincomplete_fail++; }
  }

  if($do_all || $any_fails || ($n_N > 0) || ($n_nonACGTN > 0)) { 
    $toprint = sprintf("%-*s  %-*s  %8d  %8d  %6d  %5d  %3s  %11s  %5d  %7d  %7d  ", 
                       $paccn_w+1, $paccn, $ntaccn_w, $ntaccn, $mincoord, $maxcoord, $cds_len_to_translate, $nexon, $strand2print, 
                       $incomplete, $codon_start, $n_N, $n_nonACGTN);
    #if(length($paccn)  > $paccn_w)  { $paccn_w  = length($paccn); }
    #if(length($ntaccn) > $ntaccn_w) { $ntaccn_w = length($ntaccn); }
    for(my $j = 0; $j < $last_test; $j++) { 
      $toprint .= sprintf(" %4d", $fail_A[$j]);
      if($fail_A[$j]) { 
        if($incomplete eq "no") { $complete_nfail_A[$j]++;   }
        else                    { $incomplete_nfail_A[$j]++; }
      }
    }
    $toprint .= sprintf("    %-9s", $any_fails ? "fail" : "pass");
    if($errmsg ne "") { $toprint .= " # " . $errmsg; }
    $toprint .= "\n";

    if($do_verbose) { 
      $toprint .= "TRANSLATED: $prot_translated\n"; 
      $toprint .= "FETCHED:    $prot_fetched\n"; 
    }

    if   ($ic_start == 0 && $ic_stop == 0) { push(@c_toprint_A,       $toprint); }
    elsif($ic_start == 1 && $ic_stop == 0) { push(@icstart_toprint_A, $toprint); }
    elsif($ic_start == 0 && $ic_stop == 1) { push(@icstop_toprint_A,  $toprint); }
    elsif($ic_start == 1 && $ic_stop == 1) { push(@icboth_toprint_A,  $toprint); }
  }
}
# clean up
if(-e $in_fafile . ".ssi") { unlink $in_fafile . ".ssi"; }
$sqfile->close_sqfile();

my $header_line = sprintf("%-*s  %-*s  %8s  %8s  %6s  %5s  %3s  %11s  %5s  %7s  %7s  ", 
                          $paccn_w+1, "#protein-accession", $ntaccn_w, "nt-accession", 
                          "mincoord", "maxcoord", "tr-len", "nexon", "str",
                          "incomplete?", "start", "num-Ns", "num-oth");
for(my $i = 0; $i < $last_test; $i++) { 
  $header_line .= sprintf("   T%d", ($i+1)); 
}
$header_line .= "  pass/fail\n";

# print output, first complete CDS, then incomplete
if(scalar(@c_toprint_A) > 0) { 
  print "# Complete CDS:\n";
  print $header_line;
  foreach $toprint (@c_toprint_A) { print $toprint; }
  print "#\n";
}
if(! $skip_incompletes) { 
  if(scalar(@icstart_toprint_A) > 0) { 
    print "# Incomplete CDS: (incomplete start, complete stop)\n";
    print $header_line;
    foreach $toprint (@icstart_toprint_A) { print $toprint; }
    print "#\n";
  }
  if(scalar(@icstop_toprint_A) > 0) { 
    print "# Incomplete CDS (complete start, incomplete stop):\n";
    print $header_line;
    foreach $toprint (@icstop_toprint_A) { print $toprint; }
    print "#\n";
  }
  if(scalar(@icboth_toprint_A) > 0) { 
    print "# Incomplete CDS (incomplete start and incomplete stop):\n";
    print $header_line;
    foreach $toprint (@icboth_toprint_A) { print $toprint; }
  print "#\n";
  }
}
# print summary

printf("#\n");
printf("# Summary:\n");
printf("#\n");

my $header_line2 = sprintf("%-22s  ", "# category");
for(my $i = 0; $i < $last_test; $i++) { 
  $header_line2 .= sprintf("   T%d", ($i+1)); 
}
$header_line2 .= "\n";

print $header_line2; 

printf("%-22s  ", "# num-fails-complete", "no", "-", "-", "-", "-");
for(my $j = 0; $j < $last_test; $j++) { 
  printf(" %4d", $complete_nfail_A[$j]);
}
printf("\n");

printf("%-22s  ", "# num-fails-incomplete", "yes", "-", "-", "-", "-");
for(my $j = 0; $j < $last_test; $j++) { 
  printf(" %4d", $incomplete_nfail_A[$j]);
}
printf("\n");

printf("%-22s  ", "# num-fails-all", "-", "-", "-", "-", "-");
for(my $j = 0; $j < $last_test; $j++) { 
  printf(" %4d", $complete_nfail_A[$j] + $incomplete_nfail_A[$j]);
}
printf("\n");

printf("#\n");
printf("# category    num-pass  num-fail  fract-fail\n");
printf("# complete    %8d  %8d  %10.4f\n", 
       $ncomplete - $ncomplete_fail,   
       $ncomplete_fail, 
       ($ncomplete > 0) ? ($ncomplete_fail / $ncomplete) : 0.);
printf("# incomplete  %8d  %8d  %10.4f\n", 
       $nincomplete - $nincomplete_fail,   
       $nincomplete_fail, 
       ($nincomplete > 0) ? ($nincomplete_fail / $nincomplete) : 0.);
printf("# all         %8d  %8d  %10.4f\n", 
       $nincomplete + $ncomplete - ($ncomplete_fail + $nincomplete_fail),   
       $ncomplete_fail + $nincomplete_fail, 
       (($ncomplete + $nincomplete) > 0) ? (($ncomplete_fail + $nincomplete_fail) / ($ncomplete + $nincomplete)) : 0.);


# print descriptions of each test:
printf("#\n");
printf("# A '0' in a T<n> column above indicates the sequence 'passed' the test.\n");
printf("# A '1' in a T<n> column above indicates the sequence 'failed' the test.\n");
printf("#\n");
if($do_all) { 
  printf("# Data printed for all sequences [due to -a option]\n");
}
else { 
  printf("# Only data for sequences that fail >= 1 test OR have an ambiguous nucleotide in them are printed [change to all seqs with -a]\n");
}
printf("#\n");
for(my $i = 0; $i < $last_test; $i++) { 
  printf("# Test %d (T%d): %s\n", $i+1, $i+1, $testdesc_A[$i]);
} 
printf("#\n");
printf("# Other non-obvious column meanings:\n");
printf("# \"start\":     codon_start annotation, CDS translation begins at this position of CDS (usually 1)\n");
printf("# \"tr-len\":    length of CDS to be translated (annotated length - (codon_start - 1))\n");
printf("# \"num-Ns\":    number of N nucleotides in the CDS sequence\n");
printf("# \"num-oth\":   number of other ambiguous nucleotides (non-N, non-ACGT) in the CDS sequence\n");
printf("# \"pass/fail\": \'pass\' if all tests pass (have 0s in their respective columns), else \'fail\'\n");
printf("#\n");
printf("#[ok]\n");
exit 0;




##############
# SUBROUTINES 
##############
# Subroutine: characterizeCDS()
# Args:       $cds_name_or_desc:  CDS sequence name
#             $do_compare_input:  '1' if our input is from dnaorg_compare_genomes.pl ($cds_name_or_desc is cds desc)
#                                 '0' for default input ($cds_name_or_desc is name)
# Returns:    Three values:
#             (1) $accn_coord_strand_str: five pieces of information that we will eventually print
#                 about this sequence: "<nt_accn(s)>:::<min_coord>:::<max_coord>:::<nexons>:::<strand>"
#             (2) $ic_start: '1' if sequence is incomplete at 'start'
#             (3) $ic_stop:  '1' if sequence is incomplete at 'stop'
# Dies:       If we find an incomplete character ('>' or '<') at an
#             unexpected position or we can't parse the $cds_name_or_desc
#             for some other reason.
#             
sub characterizeCDS {
  my $sub_name = "characterizeCDS()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cds_name_or_desc, $do_compare_input) = (@_);

  my $ic_start        = 0;
  my $ic_stop         = 0;
  my $expected_strand = undef;
  my $orig_cds_name_or_desc = $cds_name_or_desc;
  
  # Example $cds_name_or_desc values:
  # 
  # default ($do_compare_input == 0) examples:
  # AAH62404.1:codon_start1:BC062404.1:27:893:+:
  # AAV98535.1:codon_start1:AY680453.1:<1:>489:+:
  # AFJ92633.1:codon_start2:JQ690867.1:<1:400:+:
  #
  # $do_compare_input == 1 examples:
  # NC_004004:1059:8027:+: product:pol protein protein_id:ref|NP_658990.1| codon_start:none-annotated
  # NC_004004:1059:8027:+: product:pol protein protein_id:ref|NP_658990.1| codon_start:-
  # NC_009942:97:3552:+:NC_009942:3552:3683:+: product:truncated polyprotein protein_id:ref|YP_006485882.1| codon_start:none-annotated
  # 
  if($do_compare_input) { 
    #NC_004004:1059:8027:+: product:pol protein protein_id:ref|NP_658990.1| codon_start:none-annotated
    if($cds_name_or_desc =~ s/\s+.*$//) { 
      #NC_004004:1059:8027:+:
      ; # expected, carry on
    }
    else { die "ERROR, (-incompare) problem parsing $orig_cds_name_or_desc to determine if CDS is incomplete or not"; }
  }
  else { #default
    #AAH62404.1:codon_start1:BC062404.1:27:893:+:
    if($cds_name_or_desc =~ s/^\S+\.?\d*\:codon_start\d://) { 
      #BC062404.1:27:893:+:
      ; # expected, carry on
    }
    else { die "ERROR, (default) problem parsing $orig_cds_name_or_desc to determine if CDS is incomplete or not"; }
  }
  # remainder of cds_name is sets of 4 ':' delimited tokens
  my @cds_name_A = split(":", $cds_name_or_desc);
  my $p = 0; 
  my $exon = 1;
  # determine number of exons
  if(scalar(@cds_name_A) % 4 != 0) { die "ERROR unable to parse $orig_cds_name_or_desc while trying to determine if it is incomplete"; }
  my $nexons = scalar(@cds_name_A) / 4;
  my $nt_accn_string = "";
  my $min_coord = undef;  # min position coordinate
  my $max_coord = undef;  # max position coordinate 
  for($p = 0; $p < scalar(@cds_name_A); $p+=4) { 
    my ($nt_accn, $start, $stop, $strand) = ($cds_name_A[$p], $cds_name_A[($p+1)], $cds_name_A[($p+2)], $cds_name_A[($p+3)]);
    # only add this accession to our string of all accns if we haven't already done so.
    if($nt_accn_string !~ m/\Q$nt_accn/) { 
      if($nt_accn_string ne "") { $nt_accn_string .= ","; }
      $nt_accn_string .= $nt_accn;
    }
    if($start =~ /^\<\d+$/) { 
      if($exon != 1) { die "ERROR start incomplete at exon \#$exon (should only be possible at exon 1)\n"; }
      $ic_start = 1; 
    }
    elsif($start !~ m/^\d+$/) { 
      die "ERROR, problem parsing start coordinate $start in $orig_cds_name_or_desc while trying to determine if it is incomplete"; 
    }
    if($stop =~ /^\>\d+$/) { 
      if($exon != $nexons) { die "ERROR stop incomplete at exon \#$exon of $nexons (should only be possible at final exon)\n"; }
      $ic_stop = 1; 
    }
    elsif($stop !~ m/^\d+$/) { 
      die "ERROR, problem parsing stop coordinate $stop in $orig_cds_name_or_desc while trying to determine if it is incomplete"; 
    }
    if(! defined $expected_strand) { 
      $expected_strand = $strand; 
    }
    elsif($expected_strand ne $strand) { 
      die "ERROR, problem parsing $orig_cds_name_or_desc, multiple strands read!"; 
    }
    $exon++;

    # printf("$start..$stop\n");
    $start =~ s/\>//g;
    $start =~ s/\<//g;
    $stop  =~ s/\>//g;
    $stop  =~ s/\<//g;
    if((! defined $min_coord) || ($start < $min_coord)) { $min_coord = $start; }
    if((! defined $min_coord) || ($stop  < $min_coord)) { $min_coord = $stop;  }
    if((! defined $max_coord) || ($start > $max_coord)) { $max_coord = $start; }
    if((! defined $max_coord) || ($stop  > $max_coord)) { $max_coord = $stop;  }
  }
  if($expected_strand eq "-") { # swap $ic_start and $ic_stop
    my $tmp = $ic_start;
    $ic_start = $ic_stop;
    $ic_stop  = $tmp;
  }
  elsif($expected_strand ne "+") { 
    die "ERROR, problem parsing $orig_cds_name_or_desc, no strand read"; 
  }
  my $ret_accn_coord_nexon_strand_str = $nt_accn_string . ":::" . $min_coord . ":::" . $max_coord . ":::" . $nexons . ":::" . $expected_strand;
  return ($ret_accn_coord_nexon_strand_str, $ic_start, $ic_stop);
}

# Subroutine: translateDNA()
# Args:       $cds_seq:     CDS sequence
#             $codon_start: N=1|2|3, translation of CDS starts at position N
# Returns:    Three values:
#             1: translated protein sequence as a string
#             2: number of N     nucleotides
#             3: number of non-N nucleotides
sub translateDNA { 
  my $sub_name = "translateDNA()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cds_seq, $codon_start) = (@_);
  
  if($codon_start !~ /^[123]$/) { die "ERROR in translateDNA, invalid codon_start: $codon_start"; }

  my $length = length($cds_seq);
  my $posn = $codon_start - 1;
  my $protein = "";
  my $n_N        = 0;
  my $n_nonACGTN = 0;
  while($posn < $length) { 
    my $codon = substr($cds_seq, $posn, 3);
    $protein .= translateCodon($codon);

    # determine how many Ns and other ambig chars it has
    
    my $i = 0;
    while(($i < 3) && ($posn < $length)) { 
      my $nt = substr($cds_seq, $posn, 1);
      if   ($nt eq "N")       { $n_N++; }
      elsif($nt !~ m/[ACGT]/) { $n_nonACGTN++; }
      $posn++;
      $i++;
    }
  }

  return ($protein, $n_N, $n_nonACGTN);
}

# Subroutine: translateCodon()
# Args:       $codon: the length 3 codon
# Returns:    translated amino acid for $codon
# Reference:  http://en.wikipedia.org/wiki/Nucleic_acid_notation
sub translateCodon { 
  my $sub_name = "translateCodon()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($codon) = (@_);

  if(length($codon) > 3)  { die "ERROR codon length > 3! ($codon)"; }
  if(length($codon) == 0) { die "ERROR zero-length codon!"; }
  if(length($codon) == 1) { 
    if($codon !~ m/[ACGTWSMKRYN]/) { die "ERROR unexpected nucleotide in length 1 codon $codon"; }
    return "~"; # special character for an incomplete codon that is ambiguous
  }
  if(length($codon) == 2) { 
    if($codon !~ m/[ACGTWSMKRYN]{2}/) { die "ERROR unexpected nucleotide in length 2 codon $codon"; }
    # there's a few codons that NCBI seems to translate because the first two positions determine 
    # the AA:
    if   ($codon eq "AC") { return "T"; }
    elsif($codon eq "CC") { return "P"; }
    elsif($codon eq "CG") { return "R"; }
    elsif($codon eq "CT") { return "L"; }
    elsif($codon eq "GC") { return "A"; }
    elsif($codon eq "GG") { return "G"; }
    elsif($codon eq "GT") { return "V"; }
    elsif($codon eq "TC") { return "S"; }
    return "~"; # special character for an incomplete codon that is ambiguous
  }

  # if we get here the codon is length 3
  if($codon !~ m/[ACGTWSMKRYBDHVN]{3}/) { die "ERROR unexpected nucleotide in codon $codon"; }

  if   ($codon eq "AAA") { return "K"; } 
  elsif($codon eq "AAC") { return "N"; }
  elsif($codon eq "AAG") { return "K"; }
  elsif($codon eq "AAT") { return "N"; }
  elsif($codon eq "AAR") { return "K"; } # special
  elsif($codon eq "AAY") { return "N"; } # special

  elsif($codon eq "ACA") { return "T"; }
  elsif($codon eq "ACC") { return "T"; }
  elsif($codon eq "ACG") { return "T"; }
  elsif($codon eq "ACT") { return "T"; }
  elsif($codon eq "ACW") { return "T"; } # special 
  elsif($codon eq "ACS") { return "T"; } # special 
  elsif($codon eq "ACM") { return "T"; } # special 
  elsif($codon eq "ACK") { return "T"; } # special 
  elsif($codon eq "ACR") { return "T"; } # special 
  elsif($codon eq "ACY") { return "T"; } # special 
  elsif($codon eq "ACB") { return "T"; } # special 
  elsif($codon eq "ACD") { return "T"; } # special 
  elsif($codon eq "ACH") { return "T"; } # special 
  elsif($codon eq "ACV") { return "T"; } # special 
  elsif($codon eq "ACN") { return "T"; } # special 

  elsif($codon eq "AGA") { return "R"; }
  elsif($codon eq "AGC") { return "S"; }
  elsif($codon eq "AGG") { return "R"; }
  elsif($codon eq "AGT") { return "S"; }
  elsif($codon eq "AGR") { return "R"; } # special
  elsif($codon eq "AGY") { return "S"; } # special

  elsif($codon eq "ATA") { return "I"; }
  elsif($codon eq "ATC") { return "I"; }
  elsif($codon eq "ATG") { return "M"; }
  elsif($codon eq "ATT") { return "I"; }
  elsif($codon eq "ATW") { return "I"; } # special 
  elsif($codon eq "ATM") { return "I"; } # special 
  elsif($codon eq "ATY") { return "I"; } # special 
  elsif($codon eq "ATH") { return "I"; } # special 



  elsif($codon eq "CAA") { return "Q"; }
  elsif($codon eq "CAC") { return "H"; }
  elsif($codon eq "CAG") { return "Q"; }
  elsif($codon eq "CAT") { return "H"; }
  elsif($codon eq "CAR") { return "Q"; } # special
  elsif($codon eq "CAY") { return "H"; } # special

  elsif($codon eq "CCA") { return "P"; }
  elsif($codon eq "CCC") { return "P"; }
  elsif($codon eq "CCG") { return "P"; }
  elsif($codon eq "CCT") { return "P"; }
  elsif($codon eq "CCW") { return "P"; } # special
  elsif($codon eq "CCS") { return "P"; } # special
  elsif($codon eq "CCM") { return "P"; } # special
  elsif($codon eq "CCK") { return "P"; } # special
  elsif($codon eq "CCR") { return "P"; } # special
  elsif($codon eq "CCY") { return "P"; } # special
  elsif($codon eq "CCB") { return "P"; } # special
  elsif($codon eq "CCD") { return "P"; } # special
  elsif($codon eq "CCH") { return "P"; } # special
  elsif($codon eq "CCV") { return "P"; } # special
  elsif($codon eq "CCN") { return "P"; } # special

  elsif($codon eq "CGA") { return "R"; }
  elsif($codon eq "CGC") { return "R"; }
  elsif($codon eq "CGG") { return "R"; }
  elsif($codon eq "CGT") { return "R"; }
  elsif($codon eq "CGW") { return "R"; } # special
  elsif($codon eq "CGS") { return "R"; } # special
  elsif($codon eq "CGM") { return "R"; } # special
  elsif($codon eq "CGK") { return "R"; } # special
  elsif($codon eq "CGR") { return "R"; } # special
  elsif($codon eq "CGY") { return "R"; } # special
  elsif($codon eq "CGN") { return "R"; } # special

  elsif($codon eq "CTA") { return "L"; }
  elsif($codon eq "CTC") { return "L"; }
  elsif($codon eq "CTG") { return "L"; }
  elsif($codon eq "CTT") { return "L"; }
  elsif($codon eq "CTW") { return "L"; } # special
  elsif($codon eq "CTS") { return "L"; } # special
  elsif($codon eq "CTM") { return "L"; } # special
  elsif($codon eq "CTK") { return "L"; } # special
  elsif($codon eq "CTR") { return "L"; } # special
  elsif($codon eq "CTY") { return "L"; } # special
  elsif($codon eq "CTB") { return "L"; } # special
  elsif($codon eq "CTD") { return "L"; } # special
  elsif($codon eq "CTH") { return "L"; } # special
  elsif($codon eq "CTV") { return "L"; } # special
  elsif($codon eq "CTN") { return "L"; } # special



  elsif($codon eq "GAA") { return "E"; }
  elsif($codon eq "GAC") { return "D"; }
  elsif($codon eq "GAG") { return "E"; }
  elsif($codon eq "GAT") { return "D"; }
  elsif($codon eq "GAR") { return "E"; } # special
  elsif($codon eq "GAY") { return "D"; } # special

  elsif($codon eq "GCA") { return "A"; }
  elsif($codon eq "GCC") { return "A"; }
  elsif($codon eq "GCG") { return "A"; }
  elsif($codon eq "GCT") { return "A"; }
  elsif($codon eq "GCW") { return "A"; } # special
  elsif($codon eq "GCS") { return "A"; } # special
  elsif($codon eq "GCM") { return "A"; } # special
  elsif($codon eq "GCK") { return "A"; } # special
  elsif($codon eq "GCR") { return "A"; } # special
  elsif($codon eq "GCY") { return "A"; } # special
  elsif($codon eq "GCB") { return "A"; } # special
  elsif($codon eq "GCD") { return "A"; } # special
  elsif($codon eq "GCH") { return "A"; } # special
  elsif($codon eq "GCV") { return "A"; } # special
  elsif($codon eq "GCN") { return "A"; } # special

  elsif($codon eq "GGA") { return "G"; }
  elsif($codon eq "GGC") { return "G"; }
  elsif($codon eq "GGG") { return "G"; }
  elsif($codon eq "GGT") { return "G"; }
  elsif($codon eq "GGW") { return "G"; } # special
  elsif($codon eq "GGS") { return "G"; } # special
  elsif($codon eq "GGM") { return "G"; } # special
  elsif($codon eq "GGK") { return "G"; } # special
  elsif($codon eq "GGR") { return "G"; } # special
  elsif($codon eq "GGY") { return "G"; } # special
  elsif($codon eq "GGB") { return "G"; } # special
  elsif($codon eq "GGD") { return "G"; } # special
  elsif($codon eq "GGH") { return "G"; } # special
  elsif($codon eq "GGV") { return "G"; } # special
  elsif($codon eq "GGN") { return "G"; } # special

  elsif($codon eq "GTA") { return "V"; }
  elsif($codon eq "GTC") { return "V"; }
  elsif($codon eq "GTG") { return "V"; }
  elsif($codon eq "GTT") { return "V"; }
  elsif($codon eq "GTW") { return "V"; } # special
  elsif($codon eq "GTS") { return "V"; } # special
  elsif($codon eq "GTM") { return "V"; } # special
  elsif($codon eq "GTK") { return "V"; } # special
  elsif($codon eq "GTR") { return "V"; } # special
  elsif($codon eq "GTY") { return "V"; } # special
  elsif($codon eq "GTB") { return "V"; } # special
  elsif($codon eq "GTD") { return "V"; } # special
  elsif($codon eq "GTH") { return "V"; } # special
  elsif($codon eq "GTV") { return "V"; } # special
  elsif($codon eq "GTN") { return "V"; } # special



  elsif($codon eq "TAA") { return "*"; }
  elsif($codon eq "TAC") { return "Y"; }
  elsif($codon eq "TAG") { return "*"; }
  elsif($codon eq "TAT") { return "Y"; }
  elsif($codon eq "TAR") { return "*"; } # special
  elsif($codon eq "TAY") { return "Y"; } # special

  elsif($codon eq "TCA") { return "S"; }
  elsif($codon eq "TCC") { return "S"; }
  elsif($codon eq "TCG") { return "S"; }
  elsif($codon eq "TCT") { return "S"; }
  elsif($codon eq "TCW") { return "S"; } # special
  elsif($codon eq "TCS") { return "S"; } # special
  elsif($codon eq "TCM") { return "S"; } # special
  elsif($codon eq "TCK") { return "S"; } # special
  elsif($codon eq "TCR") { return "S"; } # special
  elsif($codon eq "TCY") { return "S"; } # special
  elsif($codon eq "TCB") { return "S"; } # special
  elsif($codon eq "TCD") { return "S"; } # special
  elsif($codon eq "TCH") { return "S"; } # special
  elsif($codon eq "TCV") { return "S"; } # special
  elsif($codon eq "TCN") { return "S"; } # special

  elsif($codon eq "TGA") { return "*"; }
  elsif($codon eq "TGC") { return "C"; }
  elsif($codon eq "TGG") { return "W"; }
  elsif($codon eq "TGT") { return "C"; }
  elsif($codon eq "TGY") { return "C"; } # special 

  elsif($codon eq "TTA") { return "L"; }
  elsif($codon eq "TTC") { return "F"; }
  elsif($codon eq "TTG") { return "L"; }
  elsif($codon eq "TTT") { return "F"; }
  elsif($codon eq "TTR") { return "L"; }
  elsif($codon eq "TTY") { return "F"; } # special

  # and the really special
  elsif($codon eq "YTG") { return "L"; }
  elsif($codon eq "YTA") { return "L"; }
  elsif($codon eq "YTR") { return "L"; }

  elsif($codon eq "MGR") { return "R"; }
  elsif($codon eq "MGA") { return "R"; }
  elsif($codon eq "MGG") { return "R"; }

  # and the really really special: three AA ambiguities B (D or N), Z (Q or E), and J (L or I)
  # I found J here: http://www.ddbj.nig.ac.jp/sub/ref2-e.html
  elsif($codon eq "RAC") { return "B"; }
  elsif($codon eq "RAT") { return "B"; }
  elsif($codon eq "RAY") { return "B"; }

  elsif($codon eq "SAA") { return "Z"; }
  elsif($codon eq "SAG") { return "Z"; }
  elsif($codon eq "SAR") { return "Z"; }

  elsif($codon eq "WTA") { return "J"; }
  elsif($codon eq "YTA") { return "J"; }
  elsif($codon eq "YTG") { return "J"; }
  elsif($codon eq "MTA") { return "J"; }
  elsif($codon eq "MTC") { return "J"; }
  elsif($codon eq "MTT") { return "J"; }
  elsif($codon eq "MTM") { return "J"; }
  elsif($codon eq "MTW") { return "J"; }
  elsif($codon eq "MTY") { return "J"; }
  elsif($codon eq "MTH") { return "J"; }

  if($do_verbose) { print "translating $codon to X\n"; }
  return "X"; 
}
