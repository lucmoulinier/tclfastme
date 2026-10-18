#!/usr/bin/env tclsh
# Smoke test: random symmetric distance matrices, every method combination.
# Run from the build directory:  tclsh generic/test.tcl

package require tcltest
namespace import tcltest::*

set dir [pwd]
source [file join $dir pkgIndex.tcl]
package require tclfastme

proc randomMatrix {nb} {
    set Lnom {}
    for {set i 0} {$i < $nb} {incr i} {
	lappend Lnom "seq-$i"
	set T($i,$i) 0.0
	for {set j 0} {$j < $i} {incr j} {
	    set T($i,$j) [expr {rand()}]
	    set T($j,$i) $T($i,$j)
	}
    }
    set Ldst {}
    for {set i 0} {$i < $nb} {incr i} {
	for {set j 0} {$j < $nb} {incr j} {
	    lappend Ldst $T($i,$j)
	}
    }
    return [list $Lnom $Ldst]
}

# Number of leaves in a Newick string
proc nbLeaves {nwk} {
    llength [regexp -all -inline {[(,][^(),:]+:} $nwk]
}

set M4 {0 2 4 5  2 0 4 5  4 4 0 3  5 5 3 0}

foreach nb {10 300 100} {
    lassign [randomMatrix $nb] Lnom Ldst
    foreach init {bal ols bionj nj} {
	foreach swap {bal ols none} {
	    test random-$nb-$init-$swap "$nb taxa, -init_tree $init -swap $swap" {
		nbLeaves [fastme $Lnom $Ldst -init_tree $init -swap $swap]
	    } $nb
	}
    }
}

test known-4 {additive 4-taxa matrix} {
    fastme {a b c d} $M4
} "(b:1.000000,(c:1.000000,d:2.000000):2.000000,a:1.000000);\n"
test known-4-bionj {additive 4-taxa matrix, BIONJ} {
    nbLeaves [fastme {a b c d} $M4 -init_tree bionj -swap none -branch_lengths none]
} 4
test one {single species} {fastme {a} {0}} "a;\n"
test two {two species} -body {nbLeaves [fastme {a b} {0 1 1 0}]} -result 1
test two-bionj {two species with bionj} -body {
    fastme {a b} {0 1 1 0} -init_tree bionj
} -match glob -result {(*)*;*}
test empty {no species} {fastme {} {}} {}
test short {matrix too short} -body {fastme {a b c d} {0 1 2}} \
    -returnCodes error -result {distance list has 3 elements, expected 4 x 4 = 16}
test long {matrix too long} -body {fastme {a b c} $M4} \
    -returnCodes error -result {distance list has 16 elements, expected 3 x 3 = 9}
test nonnum {non numeric distance} -body {
    fastme {a b c d} {0 1 2 x 1 0 1 1 2 1 0 1 x 1 1 0}
} -returnCodes error -result {expected floating-point number but got "x"}
test badopt {unknown option} -body {fastme {a b c d} $M4 -foo bal} \
    -returnCodes error -result {bad option "-foo": must be -init_tree, -branch_lengths, or -swap}
test badval {bad option value} -body {fastme {a b c d} $M4 -swap zzz} \
    -returnCodes error -result {bad value "zzz" for -swap: must be bal, ols or none}
test missingval {option without value} -body {fastme {a b c d} $M4 -swap} \
    -returnCodes error -match glob -result {wrong # args: *}
test dup {duplicate names} -body {fastme {a a b c} $M4} \
    -returnCodes error -result {duplicate species name "a"}
test names {Newick-unsafe characters are replaced} {
    lsort [regexp -all -inline {[^(),:;\n]+(?=:)} \
	       [fastme {{a b} c:d (e) f;g} $M4 -branch_lengths none -swap none]]
} {_e_ a_b c_d f_g}
test case {option values are case-insensitive and not modified} {
    set v BAL
    fastme {a b c d} $M4 -swap $v -init_tree $v
    set v
} BAL

cleanupTests
