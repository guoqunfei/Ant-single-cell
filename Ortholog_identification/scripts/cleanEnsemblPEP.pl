#!/usr/bin/perl -w
use strict;
die "Usage: <Ensembl pep>\n" unless @ARGV == 1;

if ($ARGV[0] =~ /\.gz$/) {
    open IN, "gunzip -c $ARGV[0] | ";
} else {
    open IN, $ARGV[0];
}
while (<IN>) {
    if (/^>/) {
        s/>//g;
        print ">$_";
    } else {
        print;
    }
}
close IN;
