# Generate blockfile from .gro file for propane
#
# WARNING: delete headers and box size before parsing gro file.
# Box size should be set at the very beginning of this script.
# 
# This is not a generic script: it was designed for propane specifically.
# Adapt script to other systems accordingly.
#
#

# Box size - copy bottom of .gro file
setmd box_l  4.96337   5.13917   4.52386


set listA ""
set listB  ""
# if gro is loaded, read it rather than creating a new system
if { $argc > 0} {
    set grofile [lindex $argv 0]
    puts "Reading $grofile"
    
    set fp [open $grofile r]
    set file_data [read $fp]
    close $fp

    # particle index
    set i 0
    # molecule index
    set m 0
    # process data
    set data [split $file_data "\n"]
    foreach line $data {				
	if {[lindex $line 0]!=""} {
	    if { [lindex $line 1]=="A1" || [lindex $line 1]=="A2" } { 
		set type 0
		set mass 15.035
	    } elseif { [lindex $line 1]=="B1" } { 
		set type 1
		set mass 14.027
	    } else { 
		puts "Error: unknown particle type [lindex $line 1]"
		exit 1
	    }
	    part $i pos [lindex $line 3] [lindex $line 4] [lindex $line 5] \
		v [lindex $line 6] [lindex $line 7] [lindex $line 8] \
		type $type molecule $m mass $mass	    
	    # Append lists of particle indexes
	    if { [lindex $line 1]=="A1" } {
		lappend listA $i
	    } elseif { [lindex $line 1]=="A2" } {
		lappend listA $i
		incr m
	    } elseif { [lindex $line 1]=="B1" } {
		lappend listB  $i
	    } else {
		puts "Error:unkown particle type [lindex $line 1]"
		exit 1
	    }
	    incr i
	}
    }
} else {
    puts "Error: no file loaded"
    exit 1
}


setmd time 0.0
setmd skin 0.4
setmd time_step 0.001

thermostat langevin 1.0 1.0

# tabulated interaction - *generating* the blockfile requires the tabulated
# files to exist
inter 0 0 tabulated "table_A_A.tab"
inter 1 1 tabulated "table_B_B.tab"
inter 0 1 tabulated "table_A_B.tab"
# Bonded interactions
set j 0
while { $j < [setmd n_part] } {
    # Bond A1-B
    inter $j tabulated bond "table_b1.tab"
    part $j bond $j [expr $j+1]
    incr j
    # Bond B-A2
    inter $j tabulated bond "table_b1.tab"
    part $j bond $j [expr $j+1]
    incr j
    # Bond A1-B-A2
    inter $j tabulated angle "table_a1.tab"
    part [expr $j-1] bond $j [expr $j-2] $j
    incr j
}

# Create blockfile
set out [open "|gzip -c - > conf.esp.gz" w]
blockfile $out write variable all
blockfile $out write particles [list id type molecule mass pos v]
blockfile $out write interactions
blockfile $out write bonds
blockfile $out write thermostat
blockfile $out write tclvariable [list listA listB]
close $out

exit
