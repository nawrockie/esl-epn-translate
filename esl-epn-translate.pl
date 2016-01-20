#!/usr/bin/env perl
# 
# esl-epn-translate.pl:
#
#  Translate DNA sequences into protein sequences.
#
# This script uses BioEasel's SqFile module.
# 
# EPN, Wed Nov 18 13:29:09 2015

use strict;
use Getopt::Long;
use Bio::Easel::SqFile;

my $in_fafile        = "";  # name of input fasta file

# options affecting output in normal (translate) mode
my $do_translate   = 1; # changed to '0' if one of alternative output modes is selected
my $do_reqstart    = 0; # if '1' only output sequences that start with 'ATG'
my $do_reqstop     = 0; # if '1' only output sequences that stop with valid stop codon
my $do_endatstop   = 0; # if '1' stop translating at first stop encountered, changed to 1 if -endatstop is invoked
my $do_nostop      = 0; # if '1' do not translate stop to '*', changed to 1 if -endatstop is invoked

# options for alternative output:
my $do_liststarts = 0; # if '1' DO NOT translate the sequences, instead list positions/frames of all start codons
my $do_startstop  = 0; # if '1' DO NOT translate the sequences, instead find position of first in-frame stop 
                       # codon and report that for each seq, 0 for none


&GetOptions( "reqstart"    => \$do_reqstart,
             "reqstop"     => \$do_reqstop,
             "endatstop"   => \$do_endatstop,
             "nostop"      => \$do_nostop,
             "liststarts"  => \$do_liststarts,
             "startstop"   => \$do_startstop)  || die "ERROR unknown option";

my $usage;
$usage  = "esl-epn-translate.pl [OPTIONS] <input fasta file to translate (or analyze)>\n\n";
$usage .= "\tOPTIONS THAT AFFECT TRANSLATION:\n";
$usage .= "\t\t-reqstart   : only output translated sequences that start with a valid start codon\n";
$usage .= "\t\t-reqstop    : only output translated sequences that stop  with a valid stop  codon\n";
$usage .= "\t\t-endatstop  : terminate translation at first stop codon      [default: keep going]\n";
$usage .= "\t\t-nostop     : do not print stop codons (requires -endatstop) [default: do, as '*' chars]\n";
$usage .= "\tOPTIONS FOR ALTERNATIVE OUTPUT (NO TRANSLATION PERFORMED)\n";
$usage .= "\t\t-liststarts : list start positions and frames of all start codons in each sequence <position 1..L>:<frame 0|1|2>, ...\n";
$usage .= "\t\t-startstop  : output three numbers per sequence on same line, separated by space:\n";
$usage .= "\t\t            : <d1>: '1' if first 3 nt are a valid start codon (ATG), else '0'\n";
$usage .= "\t\t            : <d2>: '1' if final 3 nt are a valid in-frame stop codon (TAA|TAG|TGA), else '0'\n";
$usage .= "\t\t            : <d3>: 1st position of first in-frame stop codon [1..seqlen], 0 if none found\n";
#$usage .= "\t\t-skipinc   : skip examination of incomplete CDS'\n";

if($do_nostop && (! $do_endatstop)) { 
  die "ERROR -nostop requires -endatstop";
}
if($do_startstop && $do_liststarts) { 
  die "ERROR -startstop and -liststarts cannot be used in combination";
}

# turn off translate mode if nec
if($do_liststarts || $do_startstop) { 
  $do_translate = 0;
}

if(scalar(@ARGV) != 1) { die $usage; }
($in_fafile) = @ARGV;

if(! -s $in_fafile) { die "ERROR, $in_fafile does not exist"; }

# open sequence file, index it, and determine the number of sequences in it
my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $in_fafile });
my $nseq = $sqfile->nseq_ssi;

# for each sequence in $in_fafile:
# 1. fetch the CDS sequence from $in_fafile
# 2. translate the CDS sequence
# 3. output the protein sequence

# for each sequence in $in_fafile:
for(my $i = 0; $i < $nseq; $i++) { 

  # 1. fetch the next sequence to translate from $in_fafile
  my ($cds_name, $cds_seq) = split(/\n/, $sqfile->fetch_consecutive_seqs(1, "", -1, undef));
  $cds_name =~ s/^\>//;
  my $cds_desc = $cds_name; 
  $cds_desc =~ s/^\S+//;
  $cds_desc =~ s/^\s+//;
  $cds_name =~ s/\s+.+$//;

  if($do_translate) { 
    # 2. translate the CDS sequence or 
    my ($prot_translated, $starts_with_start, $stops_with_stop, undef) = translateDNA($cds_seq, $do_endatstop, $do_nostop);

    # 3. output the protein sequence
    my $do_not_output = 0;
    if($do_reqstart && (! $starts_with_start)) { $do_not_output = 1; }
    if($do_reqstop  && (! $stops_with_stop))   { $do_not_output = 1; }
    if(! $do_not_output) { 
      printf(">%s%s\n$prot_translated\n", $cds_name, "-translated");
    }
  }
  if($do_liststarts) { 
    $cds_seq =~ tr/a-z/A-Z/;
    $cds_seq =~ tr/U/T/;
    my $posn = 0;
    my $cds_len = length($cds_seq);
    my $nstart = 0;
    print("$cds_name");
    while($posn < ($cds_len-2)) { 
      if(substr($cds_seq, $posn, 3) eq "ATG") { 
        if($nstart == 0) { print " ";  }
        else             { print(","); }
        printf("%d:%d", $posn+1, $posn%3);
        $nstart++;
      }
      $posn++;
    }
    print("\n");
  }
  elsif($do_startstop) { 
    my ($prot_translated, $starts_with_start, $stops_with_stop, $is_full) = translateDNA($cds_seq, 1, 0); # $do_endatstop: 1, $do_nostop: 0
    my $prot_len   = length($prot_translated);
    my $final_char = substr($prot_translated, -1, 1);
    my $first_stop;
    if($final_char eq "*") { 
      $first_stop = (($prot_len-1) * 3) + 1; 
    }
    else { 
      $first_stop = 0; 
    }
    printf("$cds_name $starts_with_start $stops_with_stop $first_stop\n");
  }
}

if(-e $in_fafile . ".ssi") { unlink $in_fafile . ".ssi"; }
exit 0;

##############
# SUBROUTINES 
##############

# Subroutine: translateDNA()
# Args:       $cds_seq:      CDS sequence
#             $do_endatstop: '1' to end at first stop codon
#             $do_nostop:    '1' to not output stop codon
# Returns:    Two values: 
#             $protein:           translated protein sequence as a string
#             $starts_with_start: '1' if first 3 nt are 'ATG'
#             $stops_with_stop:   '1' if final codon is 'TAA' or 'TAG' or 'TGA' that is in-frame of stop
#             $is_full:           '1' if translated protein sequence is full length 
#                                 (starts with start, stops with stop, length is a multiple of 3, no internal stops)
sub translateDNA { 
  my $sub_name = "translateDNA()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cds_seq, $do_endatstop, $do_nostop) = (@_);
  
  my $length  = length($cds_seq);
  if($length < 3) { 
    return ("", 0, 0, 0);
  }
  my $posn    = 0;
  my $protein = "";
  my $aa = undef;
  # determine if first 3 nt are a valid start
  my $start_codon = substr($cds_seq, 0,  3);
  my $stop_codon  = substr($cds_seq, -3, 3);
  my $start_aa = translateCodon($start_codon);
  my $stop_aa  = translateCodon($stop_codon);

#  my $n_N        = 0;
#  my $n_nonACGTN = 0;
  my $protein_len = 0;
  while(($length - $posn) >= 3) { 
    my $codon = substr($cds_seq, $posn, 3);
    $aa = translateCodon($codon);
    if($aa ne "*") { 
      $protein .= $aa;
      $protein_len++;
    }
    elsif(($aa eq "*") && (! $do_nostop)) { 
      $protein .= $aa;
    }
    if($do_endatstop && $aa eq "*") { 
      $posn = $length; # breaks while loop
    }
    $posn += 3;

    # block for 
    # determine how many Ns and other ambig chars it has
    # my $i = 0;
    # while(($i < 3) && ($posn < $length)) { 
    #  my $nt = substr($cds_seq, $posn, 1);
    #  if   ($nt eq "N")       { $n_N++; }
    #  elsif($nt !~ m/[ACGT]/) { $n_nonACGTN++; }
    #  $posn++;
    #  $i++;
    # }
  }

  my $starts_with_start = ($start_aa eq "M") ? 1 : 0;
  my $stops_with_stop  = (($length % 3 == 0) && ($stop_aa eq "*")) ? 1 : 0;

  # determine if we translated start..stop, and had valid start and stop
  my $is_full = 0;
  if($length % 3 == 0) { 
    if(($protein_len+1) == ($length / 3) &&  # +1 is for stop codon which we didn't include in protein_len
       $starts_with_start                && 
       $stops_with_stop) { 
      $is_full = 1;
    }
  }

  return ($protein, $starts_with_start, $stops_with_stop, $is_full);
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

  # if($do_verbose) { print "translating $codon to X\n"; }
  return "X"; 
}
