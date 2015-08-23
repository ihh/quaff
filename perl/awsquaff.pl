#!/usr/bin/env perl -w

use Net::Amazon::EC2;
use Net::Amazon::S3;
use Net::SSH;

use FindBin qw($Bin);
my $quaff = "$Bin/../bin/quaff";

my @awsArgs = qw(bucket instances ami type);
my @inputArgs = qw(prior params null ref read);
my @clientArgs = qw(saveprior saveparams savenull savecounts savecountswithprior maxiter mininc);
my @generalArgs = qw(maxiter mininc order gaporder threshold format log kmatch kmatchn kmatchband kmatchmb);
my @flagArgs = qw(force printall nothreshold noquals verbose fwdstrand global kmatchmax kmatchoff);

my %awsArg = map (($_ => undef), @awsArgs);
my %isInputArg = map (($_ => 1), @inputArgs);
my %isClientArg = map (($_ => 1), @clientArgs);
my %isGeneralArg = map (($_ => 1), @generalArgs);
my %isFlagArg = map (($_ => 1), @flagArgs);

my %serverCommand = ('align' => 'alignserver',
		     'count' => 'countserver',
		     'train' => 'countserver',
		     'overlap' => 'overlapserver');

my %clientCombo = map (($_ => 1),
		       "train -params",
		       "count -params",
		       "train -null",
		       "count -null");

my (@inputs, $command, @args, @serverArgs);
while (@ARGV) {
    my $arg = shift @ARGV;
    if ($arg =~ /^-(.*)/) {
	die "Please specify command before arguments\n" unless defined $command;
	
	my $param = $1;
	my $val;
	if (exists $awsArg{$param}) {
	    defined($awsArg{$param} = shift @ARGV) or die "$arg needs a value\n";

	} elsif ($isInputArg{$param}) {
	    defined($val = shift @ARGV) or die "$arg needs a value\n";
	    push @args, $arg, $val;
	    unless ($clientCombo{"$command $arg"}) {
		push @inputs, $val;
		push @serverArgs, $arg, $val;
	    }

	} elsif ($isClientArg{$param}) {
	    defined($val = shift @ARGV) or die "$arg needs a value\n";
	    push @args, $arg, $val;

	} elsif ($isGeneralArg{$param}) {
	    defined($val = shift @ARGV) or die "$arg needs a value\n";
	    push @args, $arg, $val;
	    push @serverArgs, $arg, $val;

	} elsif ($isFlagArg{$param} || $param =~ /^v\d*$/ || $param =~ /^v+$/) {
	    push @args, $arg;
	    push @serverArgs, $arg;

	} else {
	    die "Unrecognized argument: $arg\n";
	}
    } else {
	if (!defined($command)) {
	    $command = $arg;
	    unless (exists ($serverCommand{$command})) {
		die "Unrecognized command: $command\n";
	    }
	} else {
	    push @inputs, $arg;
	    push @args, $arg;
	    push @serverArgs, $arg;
	}
    }
}

defined($command) or die "Please specify a command {", join(",",sort(keys%serverCommand)), "}\n";

warn "Client: $quaff $command @args\n";
warn "Server: quaff $serverCommand{$command} @serverArgs\n";
warn "Inputs: @inputs\n";

# TODO: write me...
# Sync @inputs to $awsArg{'bucket'}
# Launch $awsArg{'instances'} of $awsArg{'type'} running $awsArg{'ami'}
# (figure out number of cores for that instance type)
# Wait for instances to enter running state, get public IP addresses
# ssh to instances, start quaff servers
# Run quaff job locally with remote servers
# Terminate instances
