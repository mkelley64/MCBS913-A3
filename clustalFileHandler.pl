#!/usr/bin/perl -w
#
# clustalFileHandler directory_name [output_directory_name]
#
#   Mark Kelley
#   Assignment 3
#   MCBS913 - Spring 2014
#   02/10/10
#
my $usageMsg = q(   Usage: clustalFileHandler directory_name [output_directory_name]

        Extract each sequence from a fastafile into a single string.
        Create reverse complement of sequence
        Find longest ORF

        Output for each sequence is <one line per frame>:
          
        sequenceId  frame  longestOrfLength  startPosition
          
        Output sent to standard output. );
          
use 5.10.0;
use warnings;
use strict;

use lib 'lib';


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
&checkUsage();


#++++++++++++ Main Routine +++++++++++++++++++++++++++

my $inputDir = $ARGV[0];
my $outputDir = $inputDir . "-mod";

if ($ARGV[1]) {
    $outputDir = $ARGV[1];
}



exit;

#+++++++++++++++++++++++++++++++++++++++++++++++
#                checkUsage
#
sub checkUsage()
{
    if ( @ARGV == 0 || $ARGV[0] eq "-h" ) {
        print STDERR "$usageMsg\n";
        exit;
    }
}


