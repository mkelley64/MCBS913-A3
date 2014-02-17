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
use clustalUtils qw( getGapsFromSequence
                     mergeGaps
                     getOUZCodesFromRange
                    );


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
&checkUsage();


#++++++++++++ Main Routine +++++++++++++++++++++++++++

my $inputDir = $ARGV[0];
my $outputDir = $inputDir . "-mod";

if ($ARGV[1]) {
    $outputDir = $ARGV[1];
}

# check if output directory exists.  If not, create it
if (-d $outputDir) {
    #do nothing
    say "$outputDir already exists";
} elsif (-e $outputDir) {
    # exists but is not a directory
    say "$outputDir already exists, but is not a directory";
    exit;
} else {
    mkdir $outputDir, 0755 or die "Failed to make directory $outputDir:  $!";
}

# read files from input directory
opendir (DIR, $inputDir) or die "Error opening directory $inputDir. $!";

while (my $file = readdir(DIR)) {
    # Ignore files beginning with a period
    next if ($file =~ m/^\./);
    
    processFile($inputDir, $file, $outputDir);
}

closedir(DIR);

exit;

#+++++++++++++++++++++++++++++++++++++++++++++++
#                processFile
#
sub processFile
{
    my $inDir = $_[0];
    my $inFile = $_[1];
    my $outDir = $_[2];
    
    open(IN, "$inDir/$inFile") or die "Could not open file '$inDir/$inFile' $!";
      
    # first line better be a sequence header
    my $header = <IN>;

    if (substr( $header, 0, 1 ) ne '>') {
        print "********* ERROR ********* Expecting header, got:\n $header";
        print " is this a fasta file??? ";
        return;
    }
    
    say "Processing $inFile...";
    
    my @sequences;
    my @gaps;
    my @gapLogData;
    my $maxSeqLength = 0;
    
    # for each sequence
    while ($header) {
      
        my $seq = ""; 
        my $inLine = <IN>;
        
        # read in all input lines of bases to create sequence string
        while ($inLine && substr($inLine, 0, 1 ) ne '>') {
            chomp($inLine);     # remove line feed
            $seq = $seq . $inLine;
            $inLine = <IN>;
        }
        
        # save sequence
        push(@sequences, $seq);
        
        # save max sequence length
        my $seqLength = length($seq);
        
        if ($seqLength > $maxSeqLength) {
            $maxSeqLength = $seqLength;
        }
         
        # get gaps and save
        my(@gapArray) = getGapsFromSequence($seq);
        push(@gaps, @gapArray);
        
        #--------------------------------------------------------
        $header = $inLine;    # last line read is either next header or null     
    }
    
    # merge gaps from all sequences
    my(@overlaps) = mergeGaps(\@gaps, $maxSeqLength);
    
    # check for case of no overlaps
    my $arrSize = @overlaps;
    
    if ($arrSize < 1) {
        return;
    }
    
    # get file prefix
    $inFile =~ m/(\w+)\./;
    my $filePrefix = $1;
    
    # for each overlap region detected
    for my $overlap (@overlaps) {
    
        my $ouzCode = "";
        
        for my $sequence (@sequences) {
            $ouzCode .= "\." . getOUZCodesFromRange($sequence, $overlap->{'start'}, $overlap->{'end'});
        }
        
        $ouzCode = substr($ouzCode, 1);
        
        my $logData = { 'head'     => $filePrefix,
                        'gapStart' => $overlap->{'start'},
                        'gapEnd'   => $overlap->{'end'},
                        'ouzCode'  => $ouzCode,
                        'revFlag'  => "F",
                        'comment'  => "ambiguous"
                      };
        
        push(@gapLogData, $logData);
    }
    
    #write results to log file
    writeLogFile(@gapLogData);
    
    
    # open(my $ofh, '>', "$outDir/$inFile") or die;
 
    
    # close $ofh;
    close(IN);
}


#+++++++++++++++++++++++++++++++++++++++++++++++
#               writeLogFile
#
sub writeLogFile
{
    my @logData = @_;
       
    my $logFilePrefix = getLogPrefix("clustalFileHandler");
    
    unless (open(OUTPUT, ">$logFilePrefix.log")) {
        print "Can not write to $logFilePrefix.log";
        exit;
    }
    
    for my $data (@logData) {
        print OUTPUT $data->{'head'}    .
                    "\t"                .
                    $data->{'gapStart'} .
                    "\t"                .
                    $data->{'gapEnd'}   .
                    "\t"                .
                    $data->{'ouzCode'}  .
                    "\t"                .
                    $data->{'revFlag'}  .
                    "\t"                .
                    $data->{'comment'};
    }
    
    close(OUTPUT);    
}

#+++++++++++++++++++++++++++++++++++++++++++++++
#            getLogFilePrefix
#
sub getLogFilePrefix
{
    my $inPrefix = $_;
    
    return $inPrefix;
}


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


