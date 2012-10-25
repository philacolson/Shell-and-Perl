#!/usr/bin/perl -w
use strict;

my $DebugLevel = 1;
my $delim = ", ";
my $space = " ";
CheckArgs ("Useage: $0 InFileName OutIndexFileName", "Y N");
PrintArgs() if ($DebugLevel > 1);

main();

################################################
sub main
{
	my $FileHandlePosnPrevious = 0;
	my $MaxRecordLength = 0;

	my $count = 0;
	open(DataFile, $ARGV[0]);
	perr("Finding Max Record Length");
	while (<DataFile>)
	{
		chomp;
		$count++;

		if (/^\>/)
		{
			#this version handles case where index is the only field on the line
			$_ =~ s/\>//g;
			my $RecLine = $_.$space.$FileHandlePosnPrevious;
			$MaxRecordLength = length($RecLine) if (length($RecLine) > $MaxRecordLength)
		}
		pe1("count: ", $count) if (!($count % 500000));
		## if ($count > 50000000)
		## {
			## pe1("emergency exit from while loop");
			## perr("emergency exit from while loop");
			## last;
		## }
		$FileHandlePosnPrevious = tell DataFile;
	} #while input
	close(DataFile);

	perr("MaxRecordLength: ", $MaxRecordLength);
	$count = 0;

	$FileHandlePosnPrevious = 0;
	open(DataFile, $ARGV[0]);
	open(OutFile, ">".$ARGV[1]);
	perr("Creating Index File");
	while (<DataFile>)
	{
		chomp;
		$count++;

		if (/^\>/)
		{
			#this version handles case where index is the only field on the line
			$_ =~ s/\>//g;
			my $RecLine = $_.$space.$FileHandlePosnPrevious;
			while (length($RecLine) < $MaxRecordLength)
			{
				$RecLine .= "~";
			}
	#		pe1($RecLine);
	#		pe1(length($RecLine));
			print OutFile $RecLine."\n";
		}
		pe1("count: ", $count) if (!($count % 500000));
		## if ($count > 50000000)
		## {
			## pe1("emergency exit from while loop");
			## perr("emergency exit from while loop");
			## last;
		## }
		$FileHandlePosnPrevious = tell DataFile;
	} #while input
	close(DataFile);
	close(OutFile);
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
