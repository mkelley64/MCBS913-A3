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
                    );
                    
our $VERSION = "1.00";
                    
#+++++++++++++++++++++++++++++++++++++++++++++++
#           getGapsFromSequence
#
sub getGapsFromSequence
{
    my($seq) = @_;
    
    my @gapArray;
    my $posIndex = 0;
    
    # match gap sequences
    while ($seq =~ m/([OUZ\-]+)/ig) {
        my $gap = {'start' => $-[0],
                   'gapSeq' => $1};
        
        push(@gapArray, $gap);
        
        $posIndex = pos($seq);
    }
    
    return (@gapArray);
}

1;
