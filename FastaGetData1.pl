#!/usr/bin/perl -w
use strict;

my $DebugLevel = 1;
my $delim = ", ";
my $space = " ";
CheckArgs ("Useage: $0 InFileName Index SNP_ID_File", "Y Y Y");
PrintArgs() if ($DebugLevel > 1);

my $InFileName = $ARGV[0];
my $Index = $ARGV[1];
my $SNP_ID_File = $ARGV[2];
if ($DebugLevel > 1)
{
	pe1("InFileName: ", $InFileName);
	pe1("Index: ", $Index);
	pe1("SNP_ID_File: ", $SNP_ID_File);
}
################################################

my $RecLen = 0;
my $RecCount = 0;

open(IndexFile, $Index);
SetFileVariables();
pe2 ("RecLen: ", $RecLen);
pe2 ("RecCount: ", $RecCount);

my $count = 0;
open(SNP_IDs, $SNP_ID_File);
while (<SNP_IDs>)
{
	chomp;
	pe2 ($_);
	$count++;
	my $ExternSNPid = $_;
	my $Offset = SearchItBinary($ExternSNPid);
	#pe2("Offset: ", $Offset);
	if ($Offset != -1){PrintData($Offset);}

	if ($count > 5000000)
	{
		pe1("emergency exit from while loop");
		perr("emergency exit from while loop");
		last;
	}
} #while input
close(SNP_IDs);
close(IndexFile);


################################################
sub PrintData
{
	my ($Offset) = @_;

	open(InFile, $InFileName);
	seek InFile, $Offset, 0;
	my $count = 0;
	while (<InFile>)
	{
		chomp;
		if (/^\>/ && $count != 0){last;}
		pe1 ($_);
		$count++;
	}
	close(InFile);
}

################################################
sub SearchItBinary
{
	my ($TheKey) = @_;
#	pe1 ("SearchItBinary TheKey: ", $TheKey);
#	my $RetVal = 99999;
	my $FileRec;
	my $StartRec = 0;
	my $EndRec = 0;
	my $MidPosn = 0;

	my $count = 0;
	seek IndexFile, 0, 0;
##	open(IndexFile, $Index);
	my $FoundKey = "";
	my $FoundOffset = 0;
	$EndRec = $RecCount;
	$MidPosn = int(($EndRec - $StartRec)/2 + $StartRec);
	while (($EndRec - $StartRec) > 1)
	{
		# pe1 ("StartRec, EndRec, MidPosn: ", $StartRec, $delim, $EndRec, $delim, $MidPosn);
		($FoundKey, $FoundOffset) = GetRecordAt($MidPosn);
		# pe1("found key, offset: ",$FoundKey, $delim, $FoundOffset);
		if ($FoundKey eq $TheKey)
		{
##			close(IndexFile);
			return($FoundOffset);
		}
		elsif($FoundKey gt $TheKey)
		{
			#pe1("FoundKey GT input TheKey");
			$EndRec = $MidPosn;
			$MidPosn = int(($EndRec - $StartRec)/2 + $StartRec);
		}
		else
		{
			#pe1("FoundKey LT input TheKey");
			$StartRec = $MidPosn;
			$MidPosn = int(($EndRec - $StartRec)/2 + $StartRec);
		}

		if ($count++ > $RecCount)
		{
			pe1("emergency exit from while loop");
			perr("emergency exit from while loop");
			last;
		}
	}
	($FoundKey, $FoundOffset) = GetRecordAt($RecCount);
	if ($FoundKey eq $TheKey)
	{
		return($FoundOffset);
##		close(IndexFile);
	}

##	close(IndexFile);
	return (-1);
}
	
################################################
sub GetRecordAt
{
	my ($RecNum) = @_;
	seek IndexFile, ($RecNum-1)*$RecLen, 0;
	my $DataOut = <IndexFile>;
	chomp $DataOut;
	#pe1 ("DataOut: ", $DataOut);
        my @fields = split(/\s/,$DataOut);
	$fields[1] =~ s/~//g;
	return ($fields[0], $fields[1]);
}

################################################
sub SetFileVariables
{
##	open(IndexFile, $Index);
	my $DataOut = <IndexFile>;
#	pe1("DataOut: ", $DataOut);
	if ($^O eq "MSWin32")
	{
		$RecLen = length($DataOut) + 1;
	}
	else
	{
		$RecLen = length($DataOut);
	}

	my $size = 0;
	if ($size = -s $Index){}
	else
	{
		pe1 ("size command failed");
		exit (0);
	}

	$RecCount = $size/$RecLen;
	if (($RecCount - int($RecCount)) > 0)
	{
		$RecCount = int($RecCount) + 1;
	}
##	close(IndexFile);
}

################# Begin print functions ##################
sub perr {
for (my $i=0;$i<=$#_;$i++){ print STDERR $_[$i]; }
print STDERR "\n";
}
sub peol { print ("\n");}
sub p1 { MyPrint (@_); }
sub p2 { if ($DebugLevel >= 2){ MyPrint (@_); } }
sub p3 { if ($DebugLevel == 3){ MyPrint (@_); } }
sub pe1 { MyPrint (@_); print ("\n"); }
sub pe2 { if ($DebugLevel >= 2){ MyPrint (@_); print ("\n"); } }
sub pe3 { if ($DebugLevel == 3){ MyPrint (@_); print ("\n"); } }
sub MyPrint { for (my $i=0;$i<=$#_;$i++){ print ($_[$i]); } }
################# End print (functions ##################
sub PrintArgs
{
	pe1 ("number of args: ", $#ARGV+1);
	pe1 ("ARGS ARE: ");
	for (my $i=0;$i<=$#ARGV;$i++)
	{
		pe1 ("arg $i: ",$ARGV[$i]);
	}

	pe1 ("Program file name: ", $0);
	perr ("\nProgram file name: ", $0);
	pe1 ("-------------------------\n");
}
################################################
sub CheckArgs
{
	my ($UsageMsg, $InputFileMask) = @_;
	my @InputParmArray = split(/\s/,$UsageMsg);
	if ($#InputParmArray-2 != $#ARGV)
	{
	 perr ("Incorrect parameters");
	 pe1 ("Incorrect parameters");
	 pe1 ($UsageMsg);
		PrintArgs();
	 exit (99);
	}
	my $InputFilesOK = 1;
	my @InputFileMaskArray = split(/\s/,$InputFileMask);
	for (my $i=0;$i<=$#InputFileMaskArray;$i++)
	{
		if ($InputFileMaskArray[$i] eq "Y")
		{
			if (!(-e $ARGV[$i]))
			{
				pe1("Input file not found, file name, usage: ",
					$ARGV[$i], $delim, $InputParmArray[$i+2]);
				$InputFilesOK = 0;
			}
		}
	}
	if (!$InputFilesOK)
	{
		PrintArgs();
	 exit (99);
	}
}
################################################
