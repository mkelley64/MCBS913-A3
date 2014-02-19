#!/usr/bin/perl -w
#
#   Mark Kelley
#   Assignment 3
#   MCBS913 - Spring 2014
#   02/16/10
#
use 5.10.0;
use warnings;
use strict;

use Test::More qw( no_plan );

use lib '../lib';

use clustalUtils qw( getGapsFromSequence
                     mergeGaps
                     getOUZCodesFromRange
                     getRevisedSequenceData
                    );
                    
# test data
my $seq0 = "MADTFGASCVPWMRDGPCIAQFMGREWDLLCIKJHMDFGIA";
my $seq1 = "MADTFGAU-----RDGPCIAQFMGREWDLLCIKJHMDFGIA";
my $seq2 = "MADTFGASCVPWURDGPCIAQFMGREWDLLCIKJHMDFGIA";
my $seq3 = "MADTFGA-----------UAQFMGREWDLLCIKJHMDFGIA";
my $seq4 = "MADTFGAU-----RDGPT----------ZLCIKJHMDFGIA";


# Test getGapsFromSequence
my $len;
my(@gapArray0) = getGapsFromSequence($seq0);
$len = @gapArray0;
is($len, 0, 'Handles sequence with no gaps');

my(@gapArray2) = getGapsFromSequence($seq1);
$len = @gapArray2;
is($len, 1, 'Handles sequence with a single gap, length > 1');
is($gapArray2[0]{'start'}, 7, 'Finds start location for single gap, length > 1');
is($gapArray2[0]{'end'}, 12, 'Finds sequence for single gap, length > 1');

my(@gapArray3) = getGapsFromSequence($seq4);
$len = @gapArray3;
is($len, 2, 'Handles sequence with multiple gaps');
is($gapArray3[0]{'start'}, 7, 'Finds first start location for multiple gaps');
is($gapArray3[0]{'end'}, 12, 'Finds first sequence for multiple gaps');
is($gapArray3[1]{'start'}, 18, 'Finds second start location for multiple gaps');
is($gapArray3[1]{'end'}, 28, 'Finds second sequence for multiple gaps');


# Test mergeGaps
my @testGaps;

my @testGaps0 = (
                {
                    'start' => 7,
                    'end'   => 13
                }
            );
push (@testGaps, @testGaps0);

my @testGaps1 = (
            );
push (@testGaps, @testGaps1);

my(@mergedGaps0) = mergeGaps(\@testGaps, 41);
$len = @mergedGaps0;
is($len, 0, 'Handles no overlaps');

my @testGaps2 = (
                {
                    'start' => 7,
                    'end'   => 18
                }
            );
push (@testGaps, @testGaps2);

my @testGaps3 = (
                {
                    'start' => 7,
                    'end'   => 13
                },
                {
                    'start' => 18,
                    'end'   => 28
                }
            );
push (@testGaps, @testGaps3);

my(@mergedGaps1) = mergeGaps(\@testGaps, 41);
$len = @mergedGaps1;
is($len, 1, 'Handles gap merging; single overlap');
is($mergedGaps1[0]{'start'}, 7, 'Finds start location for one overlap');
is($mergedGaps1[0]{'end'}, 28, 'Finds end location for one overlap');

my @testGaps4 = (
                {
                    'start' => 0,
                    'end'   => 3
                },
            );
push (@testGaps, @testGaps4);

my @testGaps5 = (
                {
                    'start' => 2,
                    'end'   => 4
                },
            );
push (@testGaps, @testGaps5);

my(@mergedGaps2) = mergeGaps(\@testGaps, 41);
$len = @mergedGaps2;
is($len, 2, 'Handles gap merging; multiple overlaps');
is($mergedGaps2[0]{'start'}, 0, 'Finds first start location for first overlap');
is($mergedGaps2[0]{'end'}, 4, 'Finds first end location for second overlap');
is($mergedGaps2[1]{'start'}, 7, 'Finds start location for second overlap');
is($mergedGaps2[1]{'end'}, 28, 'Finds end location for second overlap');


# Test getOUZCodesFromRange
my($ouzString) = getOUZCodesFromRange($seq1, 7, 13);
is($ouzString, "U", 'Handles one code at start of gap');

$ouzString = getOUZCodesFromRange($seq3, 7, 18);
is($ouzString, "U", 'Handles one code at end of gap');

$ouzString = getOUZCodesFromRange($seq4, 7, 28);
is($ouzString, "UZ", 'Handles multiple gaps in range');


# Test getRevisedSequenceData

#--------- Case 0 ---------------
# No alignments -> no revisions 
$seq0 = "ABCZ---";
$seq1 = "U---DEF";

my $seqArray = [$seq0, $seq1];
my ($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 0, 6);
my @revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "ABCZ---", 'Handles no revisions case, first sequence');
is($revSeqs[1], "U---DEF", 'Handles no revisions case, second sequence');
is($modFlag, "F", 'Handles no revisions case, modification flag');
is($ambigFlag, "", 'Handles no revisions case, ambiguous flag');

#--------- Case 1 ---------------
# Starts aligned, but not ends - move to start

# w/o revs
$seq0 = "Z------ABC";
$seq1 = "U---DEFGHI";
$seq2 = "O--------J";

$seqArray = [$seq0, $seq1, $seq2];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 0, 8);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "Z------ABC", 'Handles Case 1 w/o rev, first sequence');
is($revSeqs[1], "U---DEFGHI", 'Handles Case 1 w/o rev, second sequence');
is($revSeqs[2], "O--------J", 'Handles Case 1 w/o rev, third sequence');
is($modFlag, "F", 'Handles Case 1 w/o rev, modification flag');
is($ambigFlag, "", 'Handles Case 1 w/o rev, ambiguous flag');


# w/revs
$seq0 = "------ZABC";
$seq1 = "U---DEFGHI";
$seq2 = "O--------J";

$seqArray = [$seq0, $seq1, $seq2];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 0, 8);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "Z------ABC", 'Handles Case 1 w/rev, first sequence');
is($revSeqs[1], "U---DEFGHI", 'Handles Case 1 w/rev, second sequence');
is($revSeqs[2], "O--------J", 'Handles Case 1 w/rev, third sequence');
is($modFlag, "T", 'Handles Case 1 w/rev, modification flag');
is($ambigFlag, "", 'Handles Case 1 w/rev, ambiguous flag');


#--------- Case 2 --------------- 
# case where ends aligned, but not starts - move to end

# w/o revs
$seq0 = "CBA------Z";
$seq1 = "IHGFED---U";
$seq2 = "---------O";

$seqArray = [$seq0, $seq1, $seq2];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 0, 9);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "CBA------Z", 'Handles Case 2 w/o rev, first sequence');
is($revSeqs[1], "IHGFED---U", 'Handles Case 2 w/o rev, second sequence');
is($revSeqs[2], "---------O", 'Handles Case 2 w/o rev, third sequence');
is($modFlag, "F", 'Handles Case 2 w/o rev, modification flag');
is($ambigFlag, "", 'Handles Case 2 w/o rev, ambiguous flag');

# w/revs
$seq0 = "CBAZ------";
$seq1 = "IHGFED---U";
$seq2 = "O---------";

$seqArray = [$seq0, $seq1, $seq2];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 0, 9);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "CBA------Z", 'Handles Case 2 w/rev, first sequence');
is($revSeqs[1], "IHGFED---U", 'Handles Case 2 w/rev, second sequence');
is($revSeqs[2], "---------O", 'Handles Case 2 w/rev, third sequence');
is($modFlag, "T", 'Handles Case 2 w/rev, modification flag');
is($ambigFlag, "", 'Handles Case 2 w/rev, ambiguous flag');

#--------- Case 3A ---------------    
# case where both starts and ends aligned, no no-dash info

# w/o revs
$seq0 = "AU-----B";
$seq1 = "CO-----D";
$seq2 = "EZ-----F";
$seq3 = "GHIJKLMN";

$seqArray = [$seq0, $seq1, $seq2, $seq3];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 1, 6);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "AU-----B", 'Handles Case 3A w/o rev, first sequence');
is($revSeqs[1], "CO-----D", 'Handles Case 3A w/o rev, second sequence');
is($revSeqs[2], "EZ-----F", 'Handles Case 3A w/o rev, third sequence');
is($revSeqs[3], "GHIJKLMN", 'Handles Case 3A w/o rev, fourth sequence');
is($modFlag, "F", 'Handles Case 3A w/o rev, modification flag');
is($ambigFlag, "", 'Handles Case 3A w/o rev, ambiguous flag');

# w/revs
$seq0 = "A-----UB";
$seq1 = "CO-----D";
$seq2 = "E-----ZF";
$seq3 = "GHIJKLMN";

$seqArray = [$seq0, $seq1, $seq2, $seq3];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 1, 6);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "AU-----B", 'Handles Case 3A w/rev, first sequence');
is($revSeqs[1], "CO-----D", 'Handles Case 3A w/rev, second sequence');
is($revSeqs[2], "EZ-----F", 'Handles Case 3A w/rev, third sequence');
is($revSeqs[3], "GHIJKLMN", 'Handles Case 3A w/rev, fourth sequence');
is($modFlag, "T", 'Handles Case 3A w/rev, modification flag');
is($ambigFlag, "", 'Handles Case 3A w/rev, ambiguous flag');

#--------- Case 3B ---------------     
#case where [OUZ] in no-dashes at start - move to start

# w/o revs
$seq0 = "AU-----B";
$seq1 = "CO-----D";
$seq2 = "EZ-----F";
$seq3 = "GZIJKLMN";

$seqArray = [$seq0, $seq1, $seq2, $seq3];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 1, 6);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "AU-----B", 'Handles Case 3B w/o rev, first sequence');
is($revSeqs[1], "CO-----D", 'Handles Case 3B w/o rev, second sequence');
is($revSeqs[2], "EZ-----F", 'Handles Case 3B w/o rev, third sequence');
is($revSeqs[3], "GZIJKLMN", 'Handles Case 3B w/o rev, fourth sequence');
is($modFlag, "F", 'Handles Case 3B w/o rev, modification flag');
is($ambigFlag, "", 'Handles Case 3B w/o rev, ambiguous flag');

# w/revs
$seq0 = "A-----UB";
$seq1 = "CO-----D";
$seq2 = "E-----ZF";
$seq3 = "GZIJKLMN";

$seqArray = [$seq0, $seq1, $seq2, $seq3];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 1, 6);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "AU-----B", 'Handles Case 3B w/rev, first sequence');
is($revSeqs[1], "CO-----D", 'Handles Case 3B w/rev, second sequence');
is($revSeqs[2], "EZ-----F", 'Handles Case 3B w/rev, third sequence');
is($revSeqs[3], "GZIJKLMN", 'Handles Case 3B w/rev, fourth sequence');
is($modFlag, "T", 'Handles Case 3B w/rev, modification flag');
is($ambigFlag, "", 'Handles Case 3B w/rev, ambiguous flag');
        
#--------- Case 3C --------------- 
#case where [OUZ] in no-dashes at end - move to end

# w/o revs
$seq0 = "A-----UB";
$seq1 = "C-----OD";
$seq2 = "E-----ZF";
$seq3 = "GHIJKLUN";

$seqArray = [$seq0, $seq1, $seq2, $seq3];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 1, 6);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "A-----UB", 'Handles Case 3C w/o rev, first sequence');
is($revSeqs[1], "C-----OD", 'Handles Case 3C w/o rev, second sequence');
is($revSeqs[2], "E-----ZF", 'Handles Case 3C w/o rev, third sequence');
is($revSeqs[3], "GHIJKLUN", 'Handles Case 3C w/o rev, fourth sequence');
is($modFlag, "F", 'Handles Case 3C w/o rev, modification flag');
is($ambigFlag, "", 'Handles Case 3C w/o rev, ambiguous flag');

# w/revs
$seq0 = "AU-----B";
$seq1 = "CO-----D";
$seq2 = "E-----ZF";
$seq3 = "GHIJKLZN";

$seqArray = [$seq0, $seq1, $seq2, $seq3];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 1, 6);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "A-----UB", 'Handles Case 3C w/rev, first sequence');
is($revSeqs[1], "C-----OD", 'Handles Case 3C w/rev, second sequence');
is($revSeqs[2], "E-----ZF", 'Handles Case 3C w/rev, third sequence');
is($revSeqs[3], "GHIJKLZN", 'Handles Case 3C w/rev, fourth sequence');
is($modFlag, "T", 'Handles Case 3C w/rev, modification flag');
is($ambigFlag, "", 'Handles Case 3C w/rev, ambiguous flag');
        
#--------- Case 6 --------------- 
#case where [OUZ] in no-dashes at start AND end - move to start (ambiguous)

# w/o revs
$seq0 = "AU-----B";
$seq1 = "CO-----D";
$seq2 = "EZ-----F";
$seq3 = "GHIJKLUN";
$seq4 = "OZPQRSTV";

$seqArray = [$seq0, $seq1, $seq2, $seq3, $seq4];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 1, 6);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "AU-----B", 'Handles Case 3D w/o rev, first sequence');
is($revSeqs[1], "CO-----D", 'Handles Case 3D w/o rev, second sequence');
is($revSeqs[2], "EZ-----F", 'Handles Case 3D w/o rev, third sequence');
is($revSeqs[3], "GHIJKLUN", 'Handles Case 3D w/o rev, fourth sequence');
is($revSeqs[4], "OZPQRSTV", 'Handles Case 3D w/o rev, fourth sequence');
is($modFlag, "F", 'Handles Case 3C w/o rev, modification flag');
is($ambigFlag, "", 'Handles Case 3C w/o rev, ambiguous flag');

# w/revs
$seq0 = "AU-----B";
$seq1 = "C-----OD";
$seq2 = "E-----ZF";
$seq3 = "GHIJKLZN";
$seq4 = "OZPQRSTV";

$seqArray = [$seq0, $seq1, $seq2, $seq3, $seq4];
($revisedSeqsRef, $modFlag, $ambigFlag) = getRevisedSequenceData($seqArray, 1, 6);
@revSeqs = @{$revisedSeqsRef};
is($revSeqs[0], "AU-----B", 'Handles Case 3D w/rev, first sequence');
is($revSeqs[1], "CO-----D", 'Handles Case 3D w/rev, second sequence');
is($revSeqs[2], "EZ-----F", 'Handles Case 3D w/rev, third sequence');
is($revSeqs[3], "GHIJKLZN", 'Handles Case 3D w/rev, fourth sequence');
is($revSeqs[4], "OZPQRSTV", 'Handles Case 3D w/rev, fourth sequence');
is($modFlag, "T", 'Handles Case 3D w/rev, modification flag');
is($ambigFlag, "ambiguous", 'Handles Case 3D w/rev, ambiguous flag');

