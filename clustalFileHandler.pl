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

        Read all clustal files in a specified directory.
        Look for overlaps in sequence gaps.
        If possible, modify sequences to align gaps.
        Write revised clustal file to output directory.
        Write log file.

         );
          
use 5.10.0;
use warnings;
use strict;

use lib 'lib';
use clustalUtils qw( getGapsFromSequence
                     mergeGaps
                     getOUZCodesFromRange
                     getRevisedSequenceData
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

my @gapLogData;

while (my $file = readdir(DIR)) {
    # Ignore files beginning with a period
    next if ($file =~ m/^\./);
    
    processFile($inputDir, $file, $outputDir);
}

closedir(DIR);

#write results to log file
writeLogFile(@gapLogData);

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
        
        chomp($header);
        
        my $seqData = { 'seq'       => $seq,
                        'header'    => $header
                       };
        
        # save sequence
        push(@sequences, $seqData);
        
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
    
    close(IN);
    
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
    
    my $writeNewFileFlag = 0;
    my $newSeqsRef;
    
    # for each overlap region detected
    for my $overlap (@overlaps) {
    
        my $ouzCode = "";
        
        # get OUZcode
        for my $sequenceRecord (@sequences) {
            $ouzCode .= "\." . getOUZCodesFromRange($sequenceRecord->{'seq'}, $overlap->{'start'}, $overlap->{'end'});
        }
        
        $ouzCode = substr($ouzCode, 1);
        
        # get revision data
        my ($revisedSeqsRef, $revFlag, $ambig) = getRevisedSequenceData(\@sequences, $overlap->{'start'}, $overlap->{'end'});
        
        if ($revFlag eq "T") {
            $writeNewFileFlag = 1;
            $newSeqsRef = $revisedSeqsRef;
        }
        
        # build line of data for log file
        my $logData = { 'head'     => $filePrefix,
                        'gapStart' => $overlap->{'start'},
                        'gapEnd'   => $overlap->{'end'},
                        'ouzCode'  => $ouzCode,
                        'revFlag'  => $revFlag,
                        'comment'  => $ambig
                      };
        
        push(@gapLogData, $logData);
    }
    
    # write revised file if appropriate
    if ($writeNewFileFlag) {
        writeRevisedClustalFile($newSeqsRef, $inFile, $outputDir);
        say "Writing revised $inFile";
    }
}

#+++++++++++++++++++++++++++++++++++++++++++++++
#               writeLogFile
#
sub writeLogFile
{
    my @logData = @_;
       
    my $logFilePrefix = getLogFilePrefix("clustalFileHandler", 0);
    
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
                    $data->{'comment'}  .
                    "\n";
    }
    
    close(OUTPUT);    
}

#+++++++++++++++++++++++++++++++++++++++++++++++
#            getLogFilePrefix
#
sub getLogFilePrefix
{
    no warnings 'recursion'; # turn off recursion warnings
                           # in this block _only_
                           
    my($inPrefix, $index) = @_;
    
    my $modPrefix = $index ? "$inPrefix-$index" : $inPrefix;
      
    if (-e "$modPrefix.log") {
        if ($index > 9) {
            unlink "$modPrefix.log" or warn "Could not delete $modPrefix.log";
        
        } else {
            rename("$modPrefix.log", getLogFilePrefix($inPrefix, ++$index) . ".log") or warn "Could not rename file modPrefix.log";
        }
    }
    
    return $modPrefix;
}

#+++++++++++++++++++++++++++++++++++++++++++++++
#               writeRevisedClustalFile
#
sub writeRevisedClustalFile
{
    my ($seqRef, $outputFile, $outputDir) = @_;
    
    open(my $ofh, '>', "$outputDir/$outputFile") or die "Could not write output file $outputDir/$outputFile.  $!";
    
    my @outSeqs = @{$seqRef};
    
    for my $seqRecord (@outSeqs) {
        print $ofh "$seqRecord->{'header'}\n$seqRecord->{'seq'}\n";
    }
    
    close $ofh;
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


