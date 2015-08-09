#!/usr/bin/env perl -w

use File::Temp;

die "Usage: $0 <prog> <input> <expected>" unless @ARGV == 3;
my ($prog, $input, $expected) = @ARGV;

my $fh = File::Temp->new();
my $fname = $fh->filename;

system "$prog $input >$fname";
my $diff = `diff $fname $expected`;

if (length $diff) {
    print "`$prog $input` does not match $expected:\n";
    print `diff -y $fname $expected`;
    print "not ok\n";
} else {
    print "ok: `$prog $input` matches $expected\n";
}
