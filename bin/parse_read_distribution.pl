#!/usr/bin/env perl

#
# Purpose: change the output of read_distribution.py into field-structured data.
# See __DATA__ section for an example of such output.
#

use strict;
use warnings;

my @features = qw(Reads Tags Assigned CDS_Exons 5pUTR_Exons 3pUTR_Exons Introns TSS_up_1kb TSS_up_5kb TSS_up_10kb TES_down_1kb TES_down_5kb TES_down_10kb);
my %features;
@features{@features} = (0) x @features;


local $, = "\n";

my $sample = shift @ARGV || die "I need a sample name or --header instruction";

if ($sample eq '--header') {      # this is admittedly ugly but ...
  local $" = "\t";
  print "\t@features\n";
  exit 0;
}

while (<>) {                      # ... we do quite a bit of checking here.
  chomp;
  next if /^===/;
  next if /^Group\s+Total_bases\s+Tag_count\s+Tags/;
  my ($key, $val) = ("-") x 2;

  s/^5'/5p/;
  s/^3'/3p/;
  s/Assigned Tags/Assigned/;

  if (/^Total\s+(\w+)\s+(\d+)/) {
    ($key, $val) = ($1, $2);
  }
  elsif (/^(\w+)\s+(\d+)\s+(\d+)/) {
    ($key, $val) = ($1, $3);
  }
  else {
    die "Unexpected line $_";
  }
  die "unexpected feature [$key]" unless defined($features{$key});
  $features{$key} = $val;
}

my @values = map { $features{$_} } @features;

local $" = "\t";
print "$sample\t@values\n";


__DATA__
Total Reads                   16909
Total Tags                    22107
Total Assigned Tags           21770
=====================================================================
Group               Total_bases         Tag_count           Tags/Kb
CDS_Exons           74797099            15550               0.21
5'UTR_Exons         6838017             299                 0.04
3'UTR_Exons         30763312            2996                0.10
Introns             1415165992          2584                0.00
TSS_up_1kb          19351022            119                 0.01
TSS_up_5kb          88734076            194                 0.00
TSS_up_10kb         163641901           206                 0.00
TES_down_1kb        21035044            64                  0.00
TES_down_5kb        92316406            97                  0.00
TES_down_10kb       165985224           135                 0.00
=====================================================================
