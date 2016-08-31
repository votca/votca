#!/usr/bin/env perl
#
# Converts journal names from abbreviated or plain form
#
# (c) Copyright 2009 Jean-Olivier Irisson. GNU General Public License
#
#------------------------------------------------------------

#remove extra info from the file
system("rm -f literature.bib");
system("cat manual.bib | grep -v -e url -e doi -e issn -e month >> literature.bib");

# Get the list of journal names in a variable
open(JOURNALS,  "abbreviations.txt") or die "Can't open journalnames.txt: $!";
@journalnames = <JOURNALS>;

# Open the bib file in a variable
open(BIBFILE,  "literature.bib") or die "Can't open .bib: $!";
@file = <BIBFILE>;

# For each element of the list of journal names we extract the full name of the journal and the abbreviated name
foreach $journame (@journalnames) {
    # split based on " = "
    $journame =~ /^(.*)\s=\s(.*)$/;
    $full = $1;
    $abbrvPoint = $2;
    # removes points from the abbreviation of there are some
    $abbrvNoPoint = $abbrvPoint;
    $abbrvNoPoint =~ s/\.//g;

    # print "$full\n";
    # print "$abbrvPoint\n";
    # print "$abbrvNoPoint\n";

    # set what is replaced

    # Replace plain by abbrvPoint
    $before = $full;
    $after = $abbrvPoint;

    # # Replace plain by abbrvNoPoint
    # $before = $full;
    # $after = $abbrvNoPoint;

    # # Replace abbrvPoint by plain
    # $before = $abbrvPoint;
    # $after = $full;

    # # Replace abbrvPoint by abbrvNoPoint
    # $before = $abbrvPoint;
    # $after = $abbrvNoPoint;

    # # Replace abbrvNoPoint by plain
    # $before = $abbrvNoPoint;
    # $after = $full;

    # # Replace abbrvNoPoint by abbrvPoint
    # $before = $abbrvNoPoint;
    # $after = $abbrvPoint;

    # escape special characters
    $before =~ s/\(/\\\(/g;
    $before =~ s/\)/\\\)/g;
    $before =~ s/ICES/\{ICES\}/g;
    $before =~ s/SIAM/\{SIAM\}/g;
    $after =~ s/\(/\\\(/g;
    $after =~ s/\)/\\\)/g; 
    $after =~ s/ICES/\{ICES\}/g;
    $after =~ s/SIAM/\{SIAM\}/g;

    # print "$before\n";
    # print "$after\n";

    # replace
    foreach $line (@file){
        $line =~ s/= \{$before\}/= \{$after\}/g;
    }

}

# Print to stdout, can be redirected to a file with " > somefile.bib "
open (OUTFILE, '> literature_short.bib') or die "Can't open irisson_bib_short.bib: $!";
print OUTFILE @file;
close (OUTFILE);
