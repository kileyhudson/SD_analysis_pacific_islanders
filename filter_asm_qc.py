#!/usr/bin/env python


import argparse
import os
import subprocess
import numpy as np
import re

parser = argparse.ArgumentParser()
parser.add_argument("--szWgacGenomicSuperDupA", required = True )
parser.add_argument("--szWgacGenomicSuperDupB", required = True )
parser.add_argument("--szSampleNameA", required = True )
parser.add_argument("--szSampleNameB", required = True )
args = parser.parse_args()

assert args.szSampleNameA != args.szSampleNameB

szJustSedefBed = "just_" + args.szSampleNameA + ".bed"
szJustWgacBed  = "just_" + args.szSampleNameB + ".bed"
szInCommonBed  = args.szSampleNameA + "_vs_" + args.szSampleNameB + "_inCommon.bed"

if ( 'TMPDIR' not in os.environ ):
    
    szTmpDir = "/tmp/"

# this directory is to put temporary files in such a way as to not
# collide with other users or the same user running another sedef
if 'TMPDIR' in os.environ:
    TMPDIR = os.environ['TMPDIR']
else:
    TMPDIR = "/tmp/" + os.environ['USER'] + "/" + str( os.getpid() )
    szCommand = "mkdir -p " + TMPDIR
    print( "about to execute: " + szCommand )
    subprocess.call( szCommand, shell = True )


def makeBigBedFile( szBedFile, szBigBedFile, szBedType ):
    # sort
    szSorted = TMPDIR + "/" + szBedFile + ".sorted"
    szCommand = "sort -k1,1 -k2,2n --parallel=50 " + szBedFile + " >" + szSorted
    print( "about to execute: " + szCommand )
    subprocess.call( szCommand, shell = True )

    szCommand = "/home/hsiehph/gordo893/packages/ucsc/bedToBigBed -type=" + szBedType + " " + szSorted + " chrom.sizes " + szBigBedFile
    print( "about to execute: " + szCommand )
    subprocess.call( szCommand, shell = True )






szSedefFileA = TMPDIR + "/sedef_front_smaller.bed"
szSedefFileB = TMPDIR + "/sedef_front_larger.bed"

with open( args.szWgacGenomicSuperDupA, "r" ) as fSedef, open( szSedefFileA, "w" ) as fSedefFileA, open( szSedefFileB, "w" ) as fSedefFileB:
    n1Line = 0
    # avoid loading the entire file into memory at once
    while True:
        szLine = fSedef.readline()
        if ( szLine == "" ):
            nNumberOfLinesInSedef = n1Line
            break

        n1Line += 1

        if ( szLine.startswith( '#' )):
            continue

        aWords = re.split( r'\t|\n', szLine )

        # looks like:
        # chromosome(0)
        # start pos(1)
        # end pos(2)
        # orientation _ (underscore) or - is reverse, + is forward(5)
        # duplicated region:
        # chromosome(6)
        # start pos(7)
        # end pos(8)
        
        szChr1 = aWords[0]
        szChr2 = aWords[6]  # changed
        nStart1 = int( aWords[1] )
        nEnd1   = int( aWords[2] )
        nStart2 = int( aWords[7] ) # changed
        nEnd2   = int( aWords[8] ) # changed
        szName  = aWords[16]  # changed

        if ( szChr1 > szChr2 ):
            # swap
            ( szChr1, nStart1, nEnd1, szChr2, nStart2, nEnd2 ) = ( szChr2, nStart2, nEnd2, szChr1, nStart1, nEnd1 )
        elif( szChr1 == szChr2 and nStart1 > nStart2 ):
            # swap
            ( szChr1, nStart1, nEnd1, szChr2, nStart2, nEnd2 ) = ( szChr2, nStart2, nEnd2, szChr1, nStart1, nEnd1 )

        
        fSedefFileA.write( szChr1 + "\t" + str( nStart1 ) + "\t" + str( nEnd1 ) + "\t" + szName + "\t" + szChr2 + "\t" + str( nStart2 ) + "\t" + str( nEnd2 ) + "\t" + str( n1Line ) + "\n" )
        fSedefFileB.write( szChr2 + "\t" + str( nStart2 ) + "\t" + str( nEnd2 ) + "\t" + szName + "\t" + szChr1 + "\t" + str( nStart1 ) + "\t" + str( nEnd1 ) + "\t" + str( n1Line ) + "\n" )

#combine errors into 1 bed file 
# so much for sedef.  Now process WGAC output

"""
szWgacFileA = TMPDIR + "/wgac_front_smaller.bed"
szWgacFileB = TMPDIR + "/wgac_front_larger.bed"


with open( args.szWgacGenomicSuperDupB, "r" ) as fWgac, open( szWgacFileA, "w" ) as fWgacFileA, open( szWgacFileB, "w" ) as fWgacFileB:

    n1Line = 0
    while True:
        szLine = fWgac.readline()
        if ( szLine == "" ):
            nNumberOfLinesInWgac = n1Line
            break

        n1Line += 1
        #aWords = szLine.split('\t')
        aWords = re.split( r'\t|\n', szLine )

        # looks like:
        # chromosome(0)
        # start pos(1)
        # end pos(2)
        # orientation _ (underscore) or - is reverse, + is forward(5)
        # duplicated region:
        # chromosome(6)
        # start pos(7)
        # end pos(8)

        szChr1 = aWords[0]
        szChr2 = aWords[6]
        nStart1 = int( aWords[1] )
        nEnd1   = int( aWords[2] )
        nStart2 = int( aWords[7] )
        nEnd2   = int( aWords[8] )
        szName  = aWords[16]

        if ( szChr1 > szChr2 or ( szChr1 == szChr2 and ( nStart1 > nStart2 ) ) ):
            # swap columns
            ( szChr1, nStart1, nEnd1, szChr2, nStart2, nEnd2 ) = ( szChr2, nStart2, nEnd2, szChr1, nStart1, nEnd1 )

        fWgacFileA.write( szChr1 + "\t" + str( nStart1 ) + "\t" + str( nEnd1 ) + "\t" + szName + "\t" + szChr2 + "\t" + str( nStart2 ) + "\t" + str( nEnd2 ) + "\t" + str( n1Line ) + "\n" )
        fWgacFileB.write( szChr2 + "\t" + str( nStart2 ) + "\t" + str( nEnd2 ) + "\t" + szName + "\t" + szChr1 + "\t" + str( nStart1 ) + "\t" + str( nEnd1 ) + "\t" + str( n1Line ) + "\n" )
        
        


# now use bedtools intersect with 50% reciprocal overlap.
# change to use intersect -u 
szFirstHalfMatches = TMPDIR + "/firsthalfmatches.txt"
szLastHalfMatches = TMPDIR + "/lasthalfmatches.txt"

szCommand = "module load bedtools/2.29.2 && bedtools intersect -f 0.5 -F 0.5 -wa -wb -a " + szSedefFileA + " -b " + szWgacFileA + " >" + szFirstHalfMatches
print( "about to execute: " + szCommand )
subprocess.call( szCommand, shell = True )

szCommand = "module load bedtools/2.29.2 && bedtools intersect -f 0.5 -F 0.5 -wa -wb -a " + szSedefFileB + " -b " + szWgacFileB + " >" + szLastHalfMatches
print( "about to execute: " + szCommand )
subprocess.call( szCommand, shell = True )

# output looks like this:
# chr1(0)    4317495(1) 4318640(2) chrUn_NW_019933505v1:1-1138(3)     chrUn_NW_019933505v1(4)    1(5)       1138(6)    109(7)     chr1(8)    4317497(9) 4319729(10) data/align_both/0012/both0060028(11)        chr20(12)   28109766(13)        28111939(14)        1059(15)

# sedef fields
# chr1(0)    4317495(1) 4318640(2) chrUn_NW_019933505v1:1-1138(3)     chrUn_NW_019933505v1(4)    1(5)       1138(6)    109(7) <- sedef line #
# wgac fields
# chr1(8)    4317497(9) 4319729(10) data/align_both/0012/both0060028(11)        chr20(12)   28109766(13)        28111939(14)        1059(15) <- wgac line number

szFirstHalfMatchesLineNumbers = TMPDIR + "/first_half_matches_numbers.txt"
with open( szFirstHalfMatches, "r" ) as fFirstHalfMatches, open( szFirstHalfMatchesLineNumbers, "w" ) as fFirstHalfMatchesNumbers:
    while True:
        szLine = fFirstHalfMatches.readline()
        if ( szLine == "" ):
            break

        #aWords = szLine.split( '\t' )
        aWords = re.split( r'\t|\n', szLine )
        # sedef fields
        # chr1(0)    4317495(1) 4318640(2) chrUn_NW_019933505v1:1-1138(3)     chrUn_NW_019933505v1(4)    1(5)       1138(6)    109(7) <- sedef line #
        # wgac fields
        # chr1(8)    4317497(9) 4319729(10) data/align_both/0012/both0060028(11)        chr20(12)   28109766(13)        28111939(14)        1059(15) <- wgac line number

        szSedefLineNumber = aWords[7]
        szWgacLineNumber  = aWords[15]

        szSedefName = aWords[3]
        szWgacName  = aWords[11]

        fFirstHalfMatchesNumbers.write( szSedefLineNumber + "\t" + szWgacLineNumber + "\t" + szSedefName + "\t" + szWgacName + "\n" )

"""
szLastHalfMatchesLineNumbers  = TMPDIR + "/last_half_matches_numbers.txt"
with open( szLastHalfMatches, "r" ) as fLastHalfMatches, open( szLastHalfMatchesLineNumbers, "w" ) as fLastHalfMatchesNumbers:
    while True:
        szLine = fLastHalfMatches.readline()
        if ( szLine == "" ):
            break

        #aWords = szLine.split( '\t' )
        aWords = re.split( r'\t|\n', szLine )
        # sedef fields
        # chr1(0)    4317495(1) 4318640(2) chrUn_NW_019933505v1:1-1138(3)     chrUn_NW_019933505v1(4)    1(5)       1138(6)    109(7) <- sedef line #
        # wgac fields
        # chr1(8)    4317497(9) 4319729(10) data/align_both/0012/both0060028(11)        chr20(12)   28109766(13)        28111939(14)        1059(15) <- wgac line number

        szSedefLineNumber = aWords[7]
        szWgacLineNumber  = aWords[15]

        szSedefName = aWords[3]
        szWgacName  = aWords[11]

        fLastHalfMatchesNumbers.write( szSedefLineNumber + "\t" + szWgacLineNumber + "\t" + szSedefName + "\t" + szWgacName + "\n" )


# prepare to find pairs of line numbers that are in both files.  This
# indicates that the first locus is in sedef and wgac and the
# associated locus is also in sedef and wgac.  To prepare, sort the
# files first by sedef and then wgac line number.  Then we can read
# through both files in one pass looking for matches.

szFirstHalfMatchesLineNumbersSorted = szFirstHalfMatchesLineNumbers + ".sorted"
szLastHalfMatchesLineNumbersSorted  = szLastHalfMatchesLineNumbers  + ".sorted"

szCommand = "sort -k1,1n -k2,2n --parallel=50 " + szFirstHalfMatchesLineNumbers + " >" + szFirstHalfMatchesLineNumbersSorted
print( "about to execute: " + szCommand )
subprocess.call( szCommand, shell = True )

szCommand = "sort -k1,1n -k2,2n --parallel=50 " + szLastHalfMatchesLineNumbers + " >" + szLastHalfMatchesLineNumbersSorted
print( "about to execute: " + szCommand )
subprocess.call( szCommand, shell = True )


szLineNumbersOfMatchingSegDups = TMPDIR + "/line_numbers_of_seg_dups_in_common.txt"

with open( szFirstHalfMatchesLineNumbersSorted, "r" ) as fFirstHalfMatches, open( szLastHalfMatchesLineNumbersSorted, "r" ) as fLastHalfMatches, open( szLineNumbersOfMatchingSegDups, "w" ) as fLineNumbersOfSegDupsInCommon :

    bIncrementFirstHalf = True
    bIncrementLastHalf = True

    while True:

        if ( bIncrementFirstHalf ):
            szFirstHalfLine = fFirstHalfMatches.readline()
            if ( szFirstHalfLine == "" ):
                break

        if ( bIncrementLastHalf ):
            szLastHalfLine = fLastHalfMatches.readline()
            if ( szLastHalfLine == "" ):
                break

        #aFirstHalfWords = szFirstHalfLine.split( '\t' )
        aFirstHalfWords = re.split( r'\t|\n', szFirstHalfLine )
        #aLastHalfWords = szLastHalfLine.split( '\t' )
        aLastHalfWords = re.split( r'\t|\n', szLastHalfLine )

        nFirstSedefNumber = int( aFirstHalfWords[0] )
        nFirstWgacNumber  = int( aFirstHalfWords[1] )
        
        nLastSedefNumber  = int( aLastHalfWords[0] )
        nLastWgacNumber   = int( aLastHalfWords[1] )

        if ( nFirstSedefNumber == nLastSedefNumber and nFirstWgacNumber == nLastWgacNumber ):

            # found a match.  The first column has a match between
            # sedef line number nFirstSedefNumber and
            # nFirstWgacNumber.  The last column has a match between
            # line number nLastSedefNumber and nLastWgacNumber.
            # There is a match between both the first and last columns
            # of the sedef and wgac files iff the line numbers of the
            # first and last sedef line numbers are the same and the
            # line numbers of the first and last wgac line numbers are
            # the same.

            # record match
            fLineNumbersOfSegDupsInCommon.write( "{:d}\t{:d}\n".format( nFirstSedefNumber, nFirstWgacNumber ) )

            bIncrementFirstHalf = True
            bIncrementLastHalf = True
        else:
            if ( nFirstSedefNumber < nLastSedefNumber ):
                bIncrementFirstHalf = True
                bIncrementLastHalf = False
            elif( nFirstSedefNumber > nLastSedefNumber ):
                bIncrementFirstHalf = False
                bIncrementLastHalf = True
            elif( ( nFirstSedefNumber == nLastSedefNumber ) and ( nFirstWgacNumber < nLastWgacNumber ) ):
                bIncrementFirstHalf = True
                bIncrementLastHalf  = False
            elif( ( nFirstSedefNumber == nLastSedefNumber ) and ( nFirstWgacNumber > nLastWgacNumber ) ):
                bIncrementFirstHalf = False
                bIncrementLastHalf  = True
                
            else:
                assert False


# not using [0] since the line numbers are 1-based

aWgacLines = np.zeros( nNumberOfLinesInWgac + 1, dtype = bool)
aSedefLines = np.zeros( nNumberOfLinesInSedef + 1, dtype = bool )

with open( szLineNumbersOfMatchingSegDups, "r" ) as fLineNumbersOfSegDupsInCommon :
    
    while True:
        szLine = fLineNumbersOfSegDupsInCommon.readline()
        if ( szLine == "" ):
            break

        #aWords = szLine.split('\t')
        aWords = re.split( r'\t|\n', szLine )

        n1SedefLine = int( aWords[0] )
        n1WgacLine  = int( aWords[1] )

        aSedefLines[ n1SedefLine ] = True
        aWgacLines[ n1WgacLine ] = True

# write the sedef lines that are not matched with a wgac line

with open( args.szWgacGenomicSuperDupA, "r" ) as fWholeSedefFile, open( szJustSedefBed, "w" ) as fJustSedef:

    n1Line = 0
    while True:
        szLine = fWholeSedefFile.readline()
        if ( szLine == "" ):
            break

        n1Line += 1

        if ( szLine.startswith('#' ) ):
            continue

        if ( not aSedefLines[ n1Line ] ):
            fJustSedef.write( szLine )

with open( args.szWgacGenomicSuperDupB, "r" ) as fWholeWgacFile, open( szJustWgacBed, "w" ) as fJustWgac, open( szInCommonBed, "w" ) as fInCommon:
    
    n1Line = 0
    while True:
        szLine = fWholeWgacFile.readline()
        if ( szLine == "" ):
            break

        n1Line += 1

        if ( aWgacLines[ n1Line ] ):
            fInCommon.write( szLine )
        else:
            fJustWgac.write( szLine )

