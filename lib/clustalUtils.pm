#
#   Mark Kelley
#   Assignment 3
#   MCBS913 - Spring 2014
#   02/16/10
#
package clustalUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw( getGapsFromSequence
                     mergeGaps
                     getOUZCodesFromRange
                    );
                    
our $VERSION = "1.00";

                    
#+++++++++++++++++++++++++++++++++++++++++++++++
#           getGapsFromSequence
#
sub getGapsFromSequence
{
    my($seq) = @_;
    
    my @gapArray;
    
    # match gap sequences
    while ($seq =~ m/([OUZ]\-+|\-+[OUZ])/ig) {    
        my $gap = { 'start' => $-[0],
                    'end' => $+[0]-1
                  };
        
        push(@gapArray, $gap);
    }
    
    return (@gapArray);
}


#+++++++++++++++++++++++++++++++++++++++++++++++
#           mergeGaps
#
sub mergeGaps
{
    my($gaps_ref, $seqLength) = @_;
    
    my @gaps= @{$gaps_ref};
    
    my @db;
    
    # initialize all values to zero
    for (my $i = 0; $i < $seqLength; $i++) {
        $db[$i] = 0;
    }
    
    for my $gap (@gaps) {
        for (my $j = $gap->{'start'}; $j <= $gap->{'end'}; $j++) {
            if ($db[$j] != 2) {
                $db[$j]++;
            }
            
        }
    }
    my $joinedCounter = join('', @db);
    
    my @gapArray;
    
    # match gap sequences
    while ($joinedCounter=~ m/[12]*21*/ig) {    
        my $gap = { 'start' => $-[0],
                    'end' => $+[0]-1 
                  };
        
        push(@gapArray, $gap);
    }
    
    
    return (@gapArray);
}


#+++++++++++++++++++++++++++++++++++++++++++++++
#           getOUZCodesFromRange
#
sub getOUZCodesFromRange
{
    my($seq, $start, $end) = @_;
    
    my $OUZString = "";
    
    my $subSeq = substr( $seq, $start, $end-$start+1);
    
    while ($subSeq=~ m/([OUZ])\-+|\-([OUZ])/ig) {
        if (defined $1) {
            $OUZString .= $1;
        }
        
        if (defined $2) {
            $OUZString .= $2;
        }
    }
    
    return $OUZString;
}

1;
