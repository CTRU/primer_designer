#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (22 Oct 2014), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my @T;

my $from_column = shift || 1;
my $to_column   = 2024;

$to_column = $from_column + 1022;

while(<>) {
  chomp;
  
  my @F = split("\t", $_);


  push @T, [ $F[0], @F[$from_column..$to_column] ];
}


#print Dumper( \@T );

for(my $i = 0 ; $i < @T; $i++) {
  my $empties = 0;
  for(my $j = 0 ; $j < int( @{$T[ $i ]}); $j++) {

    if ( ! $T[ $i ][ $j ] || $T[ $i ][ $j ] eq " " || $T[ $i ][ $j ] eq "-") {
      $T[ $i ][ $j ] = "";
      $empties++;
    }
  }

  my $line = join("\t", @{$T[ $i ]});
  

  print "$line\n"  if ( $empties < int( @{ $T[ $i ]}) - 1 );
}
  
