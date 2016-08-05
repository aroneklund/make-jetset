#
#
# getFastaLengths.pl
#
# 2005/01/24
# Aron
#
# counts number of characters in each fasta record
#

my $header = '';
my $count = 0;

while( <> ) {
  chomp;
  if( /^>(.+)$/ ) {
    if( $header ne '' ) { print "$header\t$count\n"; }
    $header = $1;
    $count = 0;
  } else {
    s/\*$//; # remove trailing "stop" if necessary
    $count += length($_);
  }
}

if( $header ne '' ) { print "$header\t$count\n"; }

