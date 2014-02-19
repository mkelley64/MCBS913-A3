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
                     getRevisedSequenceData
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

#+++++++++++++++++++++++++++++++++++++++++++++++
#           getRevisedSequenceData
#
#$seq0 = "------ZABC";
#$seq1 = "U---DEFGHI";
#$seq1 = "O--------J";
sub getRevisedSequenceData
{
    my($seqs_ref, $overlapStart, $overlapEnd) = @_;
    my @inSeqs= @{$seqs_ref};

    my $modFlag = "N";
    my $ambigFlag = "";
    my $possAmbigFlag = 0;
    
    my @outSeqs;
    
    my $startOverlapMatch = 1;
    my $endOverlapMatch   = 1;
    my $noDashStartMatch  = 0;
    my $noDashEndMatch    = 0;
    
    for my $seq (@inSeqs) {
        
        my $overlapSeq = substr($seq, $overlapStart, $overlapEnd-$overlapStart + 1);
        my $doesContainDashes = (index($overlapSeq, '-') != -1);
        my $startChar = substr($seq, $overlapStart, 1);
        
        if (not $startChar =~ m/[OUZ\-]/) { 
            if ($doesContainDashes) {
                $startOverlapMatch = 0;
            }
        
        } elsif (not $doesContainDashes) {
            $noDashStartMatch = 1;
        }
        
        my $lastChar = substr($seq, $overlapEnd, 1);
        
        if (not $lastChar =~ m/[OUZ\-]/) {
            if ($doesContainDashes) {
                $endOverlapMatch = 0;
            }
            
        } elsif (not $doesContainDashes) {
            $noDashEndMatch = 1;
        }
    }
    
    my $revAction = 0;
    
    # Case 0: no alignments
    if (!$startOverlapMatch && !$endOverlapMatch) {
        return(\@inSeqs, $modFlag, $ambigFlag);
    }
      
    # Case 1: Starts aligned, but not ends - move to start
    elsif ($startOverlapMatch && !$endOverlapMatch) {
        $revAction = 1;
    }
    
    # Case 2: Ends aligned, but not starts - move to end
    elsif (!$startOverlapMatch && $endOverlapMatch) {
        $revAction = 2;
    }
    
    # Case 3: Both starts and ends aligned
    else {
    
        #case where no-dashes have no [OUZ]s - move to start
        #case where [OUZ] in no-dashes at start - move to start
        if (!$noDashEndMatch) {
            $revAction = 1;
        }
        
        #case where [OUZ] in no-dashes at end - move to end
        elsif (!$noDashStartMatch && $noDashEndMatch) {
            $revAction = 2;
        }
        
        #case where [OUZ] in no-dashes at start AND end - move to start (ambiguous)
        else {
            $possAmbigFlag = 1;
            $revAction = 1;
        }
    }
    
    # Revision actions
    
    # move to start
    if ($revAction == 1) {
        for my $seq (@inSeqs) {
            
            my $newSeq;
            my $overlapSeq = substr($seq, $overlapStart, $overlapEnd-$overlapStart + 1);
            
            if ($overlapSeq =~ s/(-+)([OUZ])/$2$1/) {
                $newSeq = substr($seq, 0, $overlapStart) .
                          $overlapSeq .
                          substr($seq, $overlapEnd + 1);  
                $modFlag = "Y";
                
                if ($possAmbigFlag) {
                    $ambigFlag = "ambiguous";
                }
            
            } else {
                $newSeq = $seq;
            }
            
            push (@outSeqs, $newSeq);
        }
    }
    
    # move to end
    elsif ($revAction == 2) {
        for my $seq (@inSeqs) {
            
            my $newSeq;
            my $overlapSeq = substr($seq, $overlapStart, $overlapEnd-$overlapStart + 1);
            
            if ($overlapSeq =~ s/([OUZ])(-+)/$2$1/) {
                $newSeq = substr($seq, 0, $overlapStart) .
                          $overlapSeq .
                          substr($seq, $overlapEnd + 1);  
                $modFlag = "Y";
            
            } else {
                $newSeq = $seq;
            }
            
            push (@outSeqs, $newSeq);
        }
    }
    
    return (\@outSeqs, $modFlag, $ambigFlag);
}


1;
