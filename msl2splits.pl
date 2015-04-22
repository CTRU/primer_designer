#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (22 Apr 2015), contact: kbr@brugger.dk

use strict;
use warnings;
use Data::Dumper;

my $msl_file = shift;
my $TRIM_LENGTH = 0;

my $seqs = readin_msl_file( $msl_file );

#print Dumper( $seqs );

my $seqs_matrix = seqs_into_matrix( $seqs );

my %kmers;
make_kmers( $seqs, 14);
make_kmers( $seqs, 16);
make_kmers( $seqs, 18);
make_kmers( $seqs, 20);
make_kmers( $seqs, 22);




#print Dumper( \%kmers );


# 
# 
# 
# Kim Brugger (22 Apr 2015)
sub make_kmers {
  my ( $seqs_hash, $kmer_length ) = @_;


  foreach my $seq_name ( keys %$seqs_hash ) {

    for(my $i =0; $i < length($$seqs_hash{ $seq_name }) - $kmer_length; $i++) {
      my $kmer = substr($$seqs_hash{ $seq_name }, $i, $kmer_length);
      next if ( $kmer =~ /-/);
      $kmers{ $i }{ $kmer_length }{ $kmer }++;
    }
  }
    



  
  
}

#print Dumper( $seqs_matrix );

my $splits_matrix = calc_splits( $seqs_matrix );

print_matrix( $splits_matrix );

# 
# 
# 
# Kim Brugger (22 Apr 2015)
sub print_matrix {
  my ( $matrix ) = @_;

#  print Dumper( \$matrix );

  my @lines;

  push @{$lines[  0 ]}, "Position";
  push @{$lines[  1 ]}, "Consensus";
  push @{$lines[  2 ]}, "A";
  push @{$lines[  3 ]}, "C";
  push @{$lines[  4 ]}, "G";
  push @{$lines[  5 ]}, "T";

  push @{$lines[  6 ]}, "14mers";
  push @{$lines[  7 ]}, "16mers";
  push @{$lines[  8 ]}, "18mers";
  push @{$lines[  9 ]}, "20mers";
  push @{$lines[ 10 ]}, "22mers";

  for( my $j =0; $j< int(@$matrix); $j++){
    push @{$lines[ 0 ]}, $j + 1;
    for( my $i =0; $i< int(@{$$matrix[0]}); $i++){

      push @{$lines[ $i + 1 ]}, $$matrix[ $j ][ $i ];
    }
    push  @{$lines[  6 ]}, int(keys %{$kmers{ $j }{ 14 }});
    push  @{$lines[  7 ]}, int(keys %{$kmers{ $j }{ 16 }});
    push  @{$lines[  8 ]}, int(keys %{$kmers{ $j }{ 18 }});
    push  @{$lines[  9 ]}, int(keys %{$kmers{ $j }{ 20 }});
    push  @{$lines[ 10 ]}, int(keys %{$kmers{ $j }{ 22 }});
  }
  
  foreach my $line ( @lines ) {
    print join("\t", @$line) . "\n";
  }  
}


# 
# 
# 
# Kim Brugger (22 Apr 2015)
sub calc_splits {
  my ( $seq_matrix) = @_;

#  print " Matrix height: " . int(@{$$seq_matrix[0]}) . "\n\n";

#  return;

  my @res;


  for( my $i =0; $i< int(@{$$seq_matrix[0]}); $i++){
#    print "$i :: ";
    my ( $A,$C,$G,$T, $total) = (0,0,0,0, 0);
    for( my $j =0; $j< int(@$seq_matrix); $j++){
      
      if ($$seq_matrix[ $j ][ $i ] eq 'A' ) {
	$A++;
      }
      if ($$seq_matrix[ $j ][ $i ] eq 'C' ) {
	$C++;
      }
      if ($$seq_matrix[ $j ][ $i ] eq 'G' ) {
	$G++;
      }
      if ($$seq_matrix[ $j ][ $i ] eq 'T' ) {
	$T++;
      }

      $total++;

#      print "$$seq_matrix[ $j ][ $i ] - ";
    }

    $res[ $i ][ 0 ] = best_base( $A,$C,$G,$T);
    $res[ $i ][ 1 ] = sprintf("%.2f", $A*100/$total);
    $res[ $i ][ 2 ] = sprintf("%.2f", $C*100/$total);
    $res[ $i ][ 3 ] = sprintf("%.2f", $G*100/$total);
    $res[ $i ][ 4 ] = sprintf("%.2f", $T*100/$total);


    

#    print "\n";
  }


#  print Dumper( \@res );

  return \@res;
  
}




# 
# 
# 
# Kim Brugger (22 Apr 2015)
sub best_base {
  my (@values )  = @_;
    
  return 'A' if ( max_value( $values[ 0 ], $values[ 1 ]) == 1 &&
		  max_value( $values[ 0 ], $values[ 2 ]) == 1 &&
		  max_value( $values[ 0 ], $values[ 3 ]) == 1);
  
  return 'C' if ( max_value( $values[ 1 ], $values[ 0 ]) == 1 &&
		  max_value( $values[ 1 ], $values[ 2 ]) == 1 &&
		  max_value( $values[ 1 ], $values[ 3 ]) == 1);

  return 'G' if ( max_value( $values[ 2 ], $values[ 0 ]) == 1 &&
		  max_value( $values[ 2 ], $values[ 1 ]) == 1 &&
		  max_value( $values[ 2 ], $values[ 3 ]) == 1);

  return 'T' if ( max_value( $values[ 3 ], $values[ 0 ]) == 1 &&
		  max_value( $values[ 3 ], $values[ 1 ]) == 1 &&
		  max_value( $values[ 3 ], $values[ 2 ]) == 1);
  


}




# 
# 
# 
# Kim Brugger (22 Apr 2015)
sub max_value {
  my ( $a, $b) = @_;
  
  return 1 if ( $a > $b);
  return 2;
  
}






# 
# 
# 
# Kim Brugger (22 Apr 2015)
sub seqs_into_matrix {
  my ( $seqs_hash ) = @_;

  my @res;

  foreach my $seq_name ( keys %$seqs_hash ) {
    

    push @res, [split("", $$seqs_hash{ $seq_name })];
  }

  return \@res;
}






# 
# 
# 
# Kim Brugger (22 Apr 2015)
sub readin_msl_file {
  my ($f ) = @_;

  my %res;
  my($name, $seq ) = (undef, undef);

  open( my $i, $f) || die "Could not open '$f': $!\n";
  while(<$i>) {
    chomp;
    
    if ( /^\>(.*)/ ) {
      if ( $name ) {
	$seq = substr( $seq, 0, $TRIM_LENGTH) if ( $TRIM_LENGTH );
	$res{ $name } = $seq;
	$seq = "";
      }

      $name = $1;
    }
    else { 
      $seq .= $_;
    }
  }

  $seq = substr( $seq, 0, $TRIM_LENGTH) if ( $TRIM_LENGTH );
  $res{ $name } = $seq  if ( $name );


  return \%res;
  
}
