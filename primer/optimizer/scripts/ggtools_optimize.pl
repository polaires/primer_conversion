#!/usr/bin/perl -w

###########################################################################################
# Overhang Optimizer
# Copyright (C) 2023 New England Biolabs, Inc.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without any implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# included NEB End User License Agreement for more details.
# 
# The above copyright notice and reference to the NEB End User License
# Agreement shall be included in all copies or substantial portions of
# the Software.
# 
# You should have received a copy of the NEB End User License Agreement
# along with this program.  If not, see <https://www.neb.com/en-us/user-license-terms-ula>.
###########################################################################################

use strict;
use Getopt::Long;

# ### MT is a better choice of the random number generator, if available.
# ### However, the default PERL RNG would work fine.
# use Math::Random::MT qw(srand rand irand);

my %SCORE = ();

my $o_size          =   -1;
my @o_olist         =   ();
my $o_mle           =   -1;
my $o_maxgc         =   -1;
my $o_maxat         =   -1;
my @exclude         =   ();
my %exclude         =   ();
my $o_iter          =  1e3;
my $o_minimize      =   "";
my $o_verbose       =    0;
my $o_batch         =    0;
my $o_eval          =   "";
my $o_init          =   "";
my $o_seqfile       =   "";
my $o_prmfile       =   "";
my @mismatch        =   ();
my %mismatch        =   ();
my $o_check_mc_move =   "";
my $o_overhang_size =    4;
my $o_ar            = 0.05;

### command-line options
GetOptions(
    "size=f"          => \$o_size,
    "olist=s"         => \@o_olist,
    "overhang_size=f" => \$o_overhang_size,
    "mle=f"           => \$o_mle,
    "maxgc=f"         => \$o_maxgc,
    "maxat=f"         => \$o_maxat,
    "batch=f"         => \$o_batch,
    "ar=f"            => \$o_ar,
    "iter=f"          => \$o_iter,
    "eval=s"          => \$o_eval,
    "exclude=s"       => \@exclude,
    "mismatch=s"      => \@mismatch,
    "minimize"        => \$o_minimize,
    "verbose+"        => \$o_verbose,
    "init=s"          => \$o_init,
    "seqfile=s"       => \$o_seqfile,
    "prmfile=s"       => \$o_prmfile,
    "check_mc_move"   => \$o_check_mc_move,
    );

### predefined list of overhangs to exclude (if any)
for my $w ( split(/,/, join(",", @exclude)) )
{
    $exclude{$w} = 1;
}

### predefined list of mismatches to ignore (if any)
for my $w ( split(/,/, join(",", @mismatch)) )
{
    my ($b1, $b2) = split(//, $w);

    $mismatch{$b1}{$b2} = 1;
    $mismatch{$b2}{$b1} = 1;
}

### command-line usage
if ( @ARGV == 0 )
{
    print "usage: $0 [options] overhang_matrix.csv\n";
    print "\n";
    print "options\n";
    print "  --size\t\tnumber of junctions (default '$o_size')\n";
    print "  --olist\t\tspecify list of overhangs (default '@o_olist')\n";
    print "  --overhang_size\toverhang size (default '$o_overhang_size')\n";
    print "  --mle\t\t\tminimum ligation efficiency (default '$o_mle')\n";
    print "  --maxgc\t\tmaximum GC content (default '$o_maxgc')\n";
    print "  --maxat\t\tmaximum AT content (default '$o_maxat')\n";
    print "  --ar\t\t\ttarget acceptance ratio (default '$o_ar')\n";
    print "  --iter\t\tnumber of iterations (default '$o_iter')\n";
    print "  --eval\t\tevaluate a single solution (default '$o_eval')\n";
    print "  --exclude\t\tcomma-separated list of overhangs to discard (default '@exclude')\n";
    print "  --mismatch\t\tcomma-separated list of mismatches to ignore in scoring (default '@mismatch')\n";
    print "  --minimize\t\tminimize ligation fidelity instead of maximizing (default '$o_minimize')\n";
    print "  --verbose\t\tincrease verbosity level (default '$o_verbose')\n";
    print "  --init\t\tcomma-separated list of overhangs to initialize MCMC (default '$o_init')\n";
    print "\n";
    exit;
}

### command-line arguments
my $ifile = shift @ARGV;

### read overhangs and associated counts
my %mat = read_overhang_matrix($ifile);

### overhangs containing mismatches to ignore (if any)
my %ignore = mismatches_to_ignore(\%mat, \%mismatch);

### load sequence (if provided)
my $sequence = "";

if ( $o_seqfile ne "" )
{
    $sequence = load_sequence($o_seqfile);
}

### decide in which mode we run
if ( $o_eval ne "" )
{
    ### evaluate a predefined set
    my @list = split(/,/, $o_eval);
    my $fidelity = ligation_fidelity(\@list);
    
    print join(",", $fidelity, @list), "\n";
}
elsif( $o_batch > 0 )
{
    ### score N randomly picked overhangs sets
    open(FH, ">", "batch_scoring.csv") || die "Can't write 'batch_scoring.csv'";
    batch_scoring($o_batch, \*FH);
    close(FH);
}
else
{
    ### find the optimal overhang set
    
    ### check for requested number of overhangs
    if ( $o_size == -1 && @o_olist == 0 )
    {
        print "[ERROR] You have to specify number of junctions using either '--size' or '--olist' options\n";
        exit;
    }

    ### get list of allowed overhangs at each junction
    my @list = ();
    my %site = ();

    if ( $o_prmfile ne "" )
    {
        ### either use param file to load all input
        my $info = process_params($o_prmfile);

        @list = @{$$info{"list"}};
        %site = %{$$info{"site"}};
    }
    else
    {
        ### or use specified overhangs from the command line
        @list = populate_overhang_list($o_size, \@o_olist, $sequence);
    }

    ### print out overhangs
    print "List of overhangs\n";
    for ( my $i = 0; $i < @list; $i++ )
    {
        print $i, "=", join(",", @{$list[$i]}), "\n";
    }

    ### find value of the scaling parameter that, on average, only a certain
    ### number of moves is accepted (defined by $o_ar; acceptance ratio)
    my @range = (find_annealing_temperature($o_ar, \@list, \%site, 1000, 0));

    # ### it is also possible to use various annealing schedules
    # my @range = annealing_range(0.95, 0.05, \@list, \%site, 1000, 0);

    ### minimization
    if ( $o_init ne "" )
    {
        ### use initialization list
        my @init = split(/,/, $o_init);

        if ( @list != @init )
        {
            print "[ERROR] Initialization list is incorrect\n";
            exit;
        }

        ### use a predefined list to init MC minimizer
        mc(\@list, \%site, \@range, $o_iter, \@init);
    }
    else
    {
        ### use best available solution so far to start MC
        if ( %SCORE )
        {
            print "Starting MC\n";

            my @init = split(/,/, $SCORE{"list"});

            print join(",",
                $SCORE{"fidelity"},
                @init,
                overhangs_to_sites(\@init, \%site)), "\n";

            mc(\@list, \%site, \@range, $o_iter, \@init);
        }
        else
        {
            mc(\@list, \%site, \@range, $o_iter, []);
        }
    }
}

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub read_overhang_matrix {
    my ($file) = @_;

    my %mat = ();
    my @col = ();

    open(FH, $file) || die "Can't open '$file'";

    while( my $line = <FH> )
    {
        chomp($line);

        my @tokens = split(/,/, $line, -1);

        my $row = shift @tokens;

        if ( @col == 0 )
        {
            @col = @tokens;
        }
        else
        {
            for ( my $i = 0; $i < @tokens; $i++ )
            {
                if ( $tokens[$i] eq "" )
                {
                    $mat{$row}{$col[$i]} = 0;
                }
                else
                {
                    $mat{$row}{$col[$i]} = $tokens[$i];
                }
            }
        }
    }

    close(FH);

    return %mat;
}

sub populate_overhang_list {
    my ($size, $olist, $sequence) = @_;

    my @temp = @$olist;
    
    if ( $size ne -1 && scalar(@temp) < $size )
    {
        push @temp, ("ALL") x ($size - scalar(@temp));
    }

    my @list = ();

    for my $spec ( @temp )
    {
        my %set = ();
        my $num = 1;
        
        if ( $spec eq "ALL" )
        {
            $spec = join(",", sort keys %mat);
        }
        elsif ( $spec =~ m/S:(\d+)-(\d+)/ )
        {
            my %overhangs = splice_sequence($sequence, $1, $2, $o_overhang_size);
            
            $spec = join(",", keys %overhangs);
        }
        
        for my $o ( split(/,/, $spec) )
        {	
            my $o_rc = reverse_complement($o);

            ### keep unique overhangs only for now (we'll add reverse complement later)
            next if ( exists $set{$o_rc} );

            ### exclude palindromic overhangs
            next if ( $o eq $o_rc );

            ### exclude specifically defined overhangs
            next if ( exists $exclude{$o} || exists $exclude{$o_rc} );

            ### exclude overhangs with low ligation efficiency
            next if ( $o_mle ne -1 && $mat{$o}{$o_rc} < $o_mle );

            ### check gc-content
            my $gc = 0;
            my $at = 0;

            for my $nt ( split(//, $o) )
            {
                if ( uc($nt) eq "G" || uc($nt) eq "C" )
                {
                    $gc++;
                }
                if ( uc($nt) eq "A" || uc($nt) eq "T" )
                {
                    $at++;
                }
            }

            ### exclude overhangs with high gc/at content
            next if ( $o_maxgc != -1 && $gc > $o_maxgc );
            next if ( $o_maxat != -1 && $at > $o_maxat );

            ### keep overhangs ordered
            $set{$o} = $num;
            $num++;
        }

        my @sorted = sort { $set{$a} <=> $set{$b} } keys %set;

        push @list, [@sorted];
    }

    return @list;
}

sub reverse_complement {
    my ($seq) = @_;
    
    my %c = ( "A" => "T",
          "C" => "G",
          "G" => "C",
          "T" => "A",
          "a" => "t",
          "c" => "g",
          "g" => "c",
          "t" => "a",
          "-" => "-" );
    
    my $newseq = join("", reverse(map { $c{$_} } split(//, $seq)));
    
    return $newseq;
}

sub mc {
    my ($list, $site, $range, $iter, $init) = @_;

    ### current solution (to avoid probing repetitve overhangs)
    my %sol = ();

    ### total number of positions
    my $size = scalar(@$list);

    ### variable positions (>1 available overhangs)
    my @posi = ();
    my $vsize = 0;

    for ( my $i = 0; $i < $size; $i++ )
    {
        if ( scalar @{$$list[$i]} > 1 )
        {
            push @posi, $i;
            $vsize++;
        }
    }

    my @result = ();
    
    ### initialize solution
    for ( my $i = 0; $i < $size; $i++ )
    {
        if ( @$init == 0 )
        {
            ### overhang list
            my @pool = @{$$list[$i]};
            
            ### pick a random overhang
            my $o = $pool[int(rand(scalar(@pool)))];

            if ( $o_check_mc_move )
            {
                ### try another overhang (TODO: DIRTY WAY TO PREVENT INFINITE LOOP)
                my $counter = 0;
                while ( exists $sol{$o}
                    || exists $sol{reverse_complement($o)}
                    || $counter > 1000 )
                {
                    $o = $pool[int(rand(scalar(@pool)))];
                    $counter++;
                }
            }

            ### store in the solution vector
            $result[$i] = $o;

            ### keep track of the current solution
            $sol{$o} = 1;
        }
        else
        {
            $result[$i] = $$init[$i];
            
            ### keep track of the current solution
            $sol{$$init[$i]} = 1;
        }
    }
    
    ### calculte fidelity of the set
    my $s0 = ligation_fidelity(\@result);

    if ( @posi == 0 )
    {
        ### nothing to optimize (no variable positions)
        # print join(",", $s0, @result), "\n";
        print join(",",
            $s0,
            @result,
            overhangs_to_sites(\@result, $site)), "\n";
        exit;
    }

    #### best achieved score
    my $best = 0;

    ### init acceptance ratio
    my $ar = 0;

    for my $exp ( @$range )
    {
        if ( $o_verbose > 1 )
        {
            print "[E = $exp]\n";
        }

        my $count = 0;
        my $accepted = 0;
        my $rejected = 0;

        while ( $count < $iter )
        {
            ### pick a random position
            my $pos = $posi[int(rand($vsize))];
            
            ### available overhangs
            my @pool = @{$$list[$pos]};

            ### pick a random overhang
            my $o = $pool[int(rand(scalar(@pool)))];

            ### try another overhang
            next if ( exists $sol{$o} );
            next if ( exists $sol{reverse_complement($o)} );
            
            ### save restore point
            my $restore = $result[$pos];

            ### update solution
            delete $sol{$restore};
            $result[$pos] = $o;
            $sol{$o} = 1;

            ### calculte fidelity of the new set
            my $s1 = ligation_fidelity(\@result);
            
            if ( ($s1 - $s0) > 0 )
            {
                ### accept move
                $s0 = $s1;
                
                ### update MC stats
                $accepted++;
                
                ### report the best result so far
                if ( ($s0 - $best) > 0 )
                {
                    $best = $s0;
                    
                    if ( $o_verbose )
                    {
                        print join(",",
                            $s1,
                            @result,
                            overhangs_to_sites(\@result, $site)), "\n";
                    }

                    if ( !%SCORE || ($s0 - $SCORE{"fidelity"}) > 0 )
                    {
                        $SCORE{"fidelity"} = $s0;
                        $SCORE{"list"} = join(",", @result);
                    }
                }
            }
            else
            {
                ### reset
                my $delta = abs($s1 - $s0);
                my $k = 1.38064852e-23;
                my $N = 6.02214085e-23;
                my $T = 2**$exp;

                my $prob = exp(-$delta / ($k * $T / $N));
                
                my $prand = rand(1.0);
                
                if ( $prand < $prob )
                {
                    ### accept move
                    $s0 = $s1;
                    
                    ### update MC stats
                    $accepted++;
                }
                else
                {
                    ### restore solution
                    delete $sol{$result[$pos]};
                    $result[$pos] = $restore;
                    $sol{$restore} = 1;

                    ### update MC stats
                    $rejected++;
                }
            }
                
            $count++;
            
            if ( $o_verbose > 2 )
            {
                printf( "%.7f,%i\r", $accepted / $count, $count );
            }
        }
        
        $ar = $accepted / $count;
    }

    return $ar;
}

sub ligation_fidelity {
    my ($list) = @_;

    ###      A1 A2   B1 B2 ...
    ###
    ### A1   11 12   11 12
    ### A2   21 22   21 22
    ###
    ### B1   11 12   .. ..
    ### B2   21 22   .. ..
    ###
    ### ...
    
    ### number of overhangs
    my $n = scalar(@$list);

    ### complimentary overhangs
    my $clist = complimentary_overhangs($list);

    my $fidelity = 1;

    for ( my $i = 0; $i < $n; $i++ )
    {
        ## 12 + 21
        my $correct += $mat{$$list[$i]}{$$clist[$i]} + $mat{$$clist[$i]}{$$list[$i]};

        my $total = 0;

        for ( my $j = 0; $j < $n; $j++ )
        {
            # ## 11 + 12 + 21 + 22
            # $total += $mat{$$list[$i]}{$$list[$j]} + $mat{$$list[$i]}{$$clist[$j]} + $mat{$$clist[$i]}{$$list[$j]} + $mat{$$clist[$i]}{$$clist[$j]};

            my $o1 = $$list[$i];
            my $o2 = $$list[$j];

            my $c1 = $$clist[$i];
            my $c2 = $$clist[$j];

            ###    o2 c2
            ### o1 11 12
            ### c1 21 22

            ### 11 + 12 + 21 + 22
            $total += $mat{$o1}{$o2} if ( $ignore{$o1}{$o2} == 0 );
            $total += $mat{$o1}{$c2} if ( $ignore{$o1}{$c2} == 0 );
            $total += $mat{$c1}{$o2} if ( $ignore{$c1}{$o2} == 0 );
            $total += $mat{$c1}{$c2} if ( $ignore{$c1}{$c2} == 0 );
        }

        if ( $total > 0 )
        {
            $fidelity *= $correct / $total;
        }
    }

    if ( $o_minimize ne "" )
    {
        return sprintf("%.15f", 1 - $fidelity);
    }

    return sprintf("%.15f", $fidelity);
}

sub annealing_range {
    my ($max, $min, $list, $site, $iter, $exp0) = @_;

    my @range = ();

    my $max_e = find_annealing_temperature(0.95, $list, $site, $iter, $exp0);
    my $min_e = find_annealing_temperature(0.05, $list, $site, $iter, $exp0);

    for ( my $i = $max_e; $i >= $min_e; $i-- )
    {
        push @range, $i;
    }

    return @range;
}

sub find_annealing_temperature {
    my ($target_ar, $list, $site, $iter, $exp0) = @_;

    print "\n";
    print "Target acceptance ratio : $target_ar\n";
    
    my %data = ();
    
    my $ar = mc($list, $site, [$exp0], $iter, []);
    $data{$exp0} = abs($target_ar - $ar);

    print "Initial acceptance ratio : $ar\n";

    if ( $ar > $target_ar )
    {
        my $exp = $exp0;

        my $counter = 0;
        
        print "Lowering E...\n";
        do {
            $exp--;
            print "  Lowering E to $exp...\n";
            $ar = mc($list, $site, [$exp], $iter, []);
            $data{$exp} = abs($target_ar - $ar);
            print "  Acceptance ratio : $ar\n";
            
            $counter++;
        }
        while ( $ar > $target_ar && $counter < 100 );
    }
    elsif ( $ar < $target_ar )
    {
        my $exp = $exp0;

        my $counter = 0;
        
        print "Increasing E...\n";
        do {
            $exp++;
            print "  Raising E to $exp...\n";
            $ar = mc($list, $site, [$exp], $iter, []);
            $data{$exp} = abs($target_ar - $ar);
            print "  Acceptance ratio : $ar\n";
            
            $counter++;
        }
        while ( $ar < $target_ar && $counter < 100 );
    }

    ### find exponent that gives acceptance ratio closest to the target
    my @sorted = sort { $data{$a} <=> $data{$b} } keys %data;
    
    return $sorted[0];
}

sub complimentary_overhangs {
    my ($list) = @_;

    my @clist = ();

    for ( my $i = 0; $i < scalar(@$list); $i++ )
    {
        push @clist, reverse_complement($$list[$i]);
    }

    return \@clist;
}

sub min {
    my ($arr) = @_;

    my $min = $$arr[0];

    for my $v ( @$arr )
    {
        if ( $min > $v )
        {
            $min = $v;
        }
    }
    
    return $min;
}

sub max {
    my ($arr) = @_;

    my $max = $$arr[0];

    for my $v ( @$arr )
    {
        if ( $max < $v )
        {
            $max = $v;
        }
    }
    
    return $max;
}

sub sum {
    my ($arr) = @_;

    my $sum = 0;

    for my $v ( @$arr )
    {
        $sum += $v;
    }
    
    return $sum;
}

sub load_sequence {
    my ($file) = @_;

    my $sequence = "";
    
    open(FA, $file) || die "Can't open '$file'";
    
    while ( my $line = <FA> )
    {
        ### remove end of line
        chomp($line);

        ### combine everything into a single string
        if ( substr($line, 0, 1) ne ">" )
        {
            $sequence .= $line;
        }
    }
    
    ### remove any white spaces
    $sequence =~ s/ //g;

    close(FA);
    
    return uc($sequence);
}

sub splice_sequence {
    my ($sequence, $region_start, $region_end, $overhang_size) = @_;

    my %overhangs = ();

    for( my $i = $region_start - 1; $i + $overhang_size <= $region_end; $i++ )
    {
        my $oh = substr($sequence, $i, $overhang_size);

        $overhangs{$oh} = 1;
    }

    return %overhangs;
}

sub process_sequence {
    my ($file, $nj) = @_;

    my @list = ();
    my $sequence = load_sequence($file);
    my %location = ();

    print length($sequence), "\n";
    print $sequence, "\n";

    return {"list" => \@list, "location" => \%location};
}

sub process_params {
    my ($file) = @_;

    my %params = ();
    
    open(PRM, $file) || die "Can't open '$file'";

    while ( my $line = <PRM> )
    {
        chomp($line);

        my ($name, $value) = split(/=/, $line);

        if ( $name eq "overhang" )
        {
            push @{$params{$name}}, [split(/,/, $value)];
        }
        else
        {
            $params{$name} = $value;
        }
    }

    close(PRM);

    my %data = ();

    for my $entry ( @{$params{"overhang"}} )
    {
        my ($id, $start, $fixed, $oh) = @$entry;

        if ( $fixed )
        {
            $data{$id}{$oh} = $start;
        }
        else
        {
            for ( my $pos = $start - 25; $pos < $start + 25; $pos++ )
            {
                my $oh = substr($params{"sequence"}, $pos, $o_overhang_size);

                # print $oh, " ", $pos, " ", abs($start - $pos), "\n";

                if ( ! exists $data{$id}{$oh} )
                {
                    $data{$id}{$oh} = $pos;
                }
                else
                {
                    if ( abs($start - $pos) < abs($start - $data{$id}{$oh}) )
                    {
                        $data{$id}{$oh} = $pos;
                    }
                }
            }
        }
    }

    my @list = (); 
    my %site = ();

    for my $id ( sort { $a <=> $b } keys %data )
    {
        for my $oh ( keys %{$data{$id}} )
        {
            push @{$list[$id]}, $oh;

            $site{$id}{$oh} = $data{$id}{$oh};
        }
    }

    return {"list" => \@list, "site" => \%site};
}

sub overhangs_to_sites {
    my ($result, $site) = @_;
    
    my @values = ();
    
    for ( my $i = 0; $i < @$result; $i++ )
    {
        if ( exists $$site{$i}{$$result[$i]} )
        {
            push @values, $$site{$i}{$$result[$i]};
        }
        else
        {
            push @values, "";
        }
    }

    return @values;
}

sub mismatches_to_ignore {
    my ($matrix, $mismatch) = @_;

    my %ignore = ();
    
    for my $o1 ( keys %$matrix )
    {
        for my $o2 ( keys %{$$matrix{$o1}} )
        {
            ### init
            $ignore{$o1}{$o2} = 0;
            
            my $o2_r = scalar reverse($o2);

            for ( my $i = 0; $i < length($o1); $i++ )
            {
                my $b1 = substr($o1, $i, 1);
                my $b2 = substr($o2_r, $i, 1);

                if ( exists $$mismatch{$b1}{$b2} || exists $$mismatch{$b2}{$b1} )
                {
                    ### update
                    $ignore{$o1}{$o2} = 1;
                }
            }
        }
    }
    
    return %ignore;
}

sub batch_scoring {
    my ($num_sol, $fh) = @_;
    
    my @list = populate_overhang_list($o_size, \@o_olist, '');

    for ( my $k = 0; $k < $num_sol; $k++ )
    {
        my %sol = ();
        my @result = ();
        
        ### initialize solution
        for ( my $i = 0; $i < @list; $i++ )
        {
            ### overhang list
            my @pool = @{$list[$i]};
            
            ### pick a random overhang
            my $o = $pool[int(rand(scalar(@pool)))];
            
            if ( $o_check_mc_move )
            {
                ### try another overhang (TODO: DIRTY WAY TO PREVENT INFINITE LOOP)
                my $counter = 0;
                while ( exists $sol{$o}
                    || exists $sol{reverse_complement($o)}
                    || $counter > 1000 )
                {
                    $o = $pool[int(rand(scalar(@pool)))];
                    $counter++;
                }
            }

            ### store in the solution vector
            $result[$i] = $o;
            
            ### keep track of the current solution
            $sol{$o} = 1;
        }
        
        my $fidelity = ligation_fidelity(\@result);
        
        print $fh join(",", $fidelity, @result), "\n";
    }
}
