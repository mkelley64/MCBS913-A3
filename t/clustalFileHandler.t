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
                    );

# Test getGapsFromSequence
my $len;
my(@gapArray0) = getGapsFromSequence("MADTFGASCVPWMRDGPCIAQFMGREWDLLCIKJHMDFGIA");
$len = @gapArray0;
is($len, 0, 'Handles sequence with no gaps');

my(@gapArray2) = getGapsFromSequence("MADTFGAU-----RDGPCIAQFMGREWDLLCIKJHMDFGIA");
$len = @gapArray2;
is($len, 1, 'Handles sequence with a single gap, length > 1');
is($gapArray2[0]{'start'}, 7, 'Finds start location for single gap, length > 1');
is($gapArray2[0]{'end'}, 12, 'Finds sequence for single gap, length > 1');

my(@gapArray3) = getGapsFromSequence("MADTFGAU-----RDGPT----------ZLCIKJHMDFGIA");
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
my($ouzString) = getOUZCodesFromRange("MADTFGAU-----RDGPCIAQFMGREWDLLCIKJHMDFGIA", 7, 13);
is($ouzString, "U", 'Handles one code at start of gap');

$ouzString = getOUZCodesFromRange("MADTFGA-----------UAQFMGREWDLLCIKJHMDFGIA", 7, 18);
is($ouzString, "U", 'Handles one code at end of gap');

$ouzString = getOUZCodesFromRange("MADTFGAU-----RDGPT----------ZLCIKJHMDFGIA", 7, 28);
is($ouzString, "UZ", 'Handles multiple gaps in range');


