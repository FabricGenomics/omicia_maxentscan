use List::Util qw[min max];
use strict;

my $datapath = $ARGV[0];

my %me2x5 = &makescorematrix($datapath.'/me2x5');
my %seq = &makesequencematrix($datapath.'/splicemodels/splice5sequences');
my $seq_length = 9;
my $varname = "";

my %bgd;
$bgd{'A'} = 0.27;
$bgd{'C'} = 0.23;
$bgd{'G'} = 0.23;
$bgd{'T'} = 0.27;


while(<STDIN>) {
    chomp;
    if (/^\s*$/) { #discard blank lines;
        next;
    } 
    elsif (/^>/) { #get header;
        $_ =~ s/\cM//g; #gets rid of carriage return
        $varname = substr $_, 1;
        next;
    }
    else {
        $_ =~ s/\cM//g; #gets rid of carriage return
        my $str = uc($_);
        my $max_score = -999;
        my $max_iteration = (length $str) - $seq_length;
        if ($max_iteration < 0) {
            die "Sequencelength too short for ".$varname;
        }
        for (my $i = 0; $i <= $max_iteration; $i++ ) {
            my $sequence = substr $str, $i, $seq_length;
            my $score = &log2(&scoreconsensus($sequence)*$me2x5{$seq{&getrest($sequence)}});
            # print $i." ".$sequence." ".$score."\n";
            $max_score = max($score, $max_score);
        }
        print sprintf("%s\t%.4f", $varname, $max_score)."\n";
    }
}

  
sub makesequencematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$_} = $n;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}
sub makescorematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$n} = $_;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

sub getrest{
  my $seq = shift;
  my @seqa = split(//,uc($seq));
  return $seqa[0].$seqa[1].$seqa[2].$seqa[5].$seqa[6].$seqa[7].$seqa[8];
}
sub scoreconsensus{
  my $seq = shift;
  my @seqa = split(//,uc($seq));
  my %bgd; 
  $bgd{'A'} = 0.27; 
  $bgd{'C'} = 0.23; 
  $bgd{'G'} = 0.23; 
  $bgd{'T'} = 0.27;  
  my %cons1;
  $cons1{'A'} = 0.004;
  $cons1{'C'} = 0.0032;
  $cons1{'G'} = 0.9896;
  $cons1{'T'} = 0.0032;
  my %cons2;
  $cons2{'A'} = 0.0034; 
  $cons2{'C'} = 0.0039; 
  $cons2{'G'} = 0.0042; 
  $cons2{'T'} = 0.9884;
  my $addscore = $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]}); 
  return $addscore;
}

sub log2{
      my ($val) = @_;
    return log($val)/log(2);
}
