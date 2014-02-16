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
                    );

# Test getGapsFromSequence
my $len;
my(@gapArray0) = getGapsFromSequence("MADTFGASCVPWMRDGPCIAQFMGREWDLLCIKJHMDFGIA");
$len = @gapArray0;
is($len, 0, 'Handles sequence with no gaps');

my(@gapArray1) = getGapsFromSequence("MADTFGASCVPWURDGPCIAQFMGREWDLLCIKJHMDFGIA");
$len = @gapArray1;
is($len, 1, 'Handles sequence with a single gap, length==1');
is($gapArray1[0]{'start'}, 12, 'Finds start location for single gap, length==1');
is($gapArray1[0]{'gapSeq'}, "U", 'Finds sequence for single gap, length==1');

my(@gapArray2) = getGapsFromSequence("MADTFGAU-----RDGPCIAQFMGREWDLLCIKJHMDFGIA");
$len = @gapArray2;
is($len, 1, 'Handles sequence with a single gap, length > 1');
is($gapArray2[0]{'start'}, 7, 'Finds start location for single gap, length > 1');
is($gapArray2[0]{'gapSeq'}, "U-----", 'Finds sequence for single gap, length > 1');