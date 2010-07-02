# Generate blockfile from .gro file for methanol
#
# WARNING: delete headers and box size before parsing gro file.
# Box size should be set at the very beginning of this script.
# 
# This is not a generic script: it was designed for methanol specifically.
# Adapt script to other systems accordingly.
#
# Usage: Espresso_bin generate.tcl conf.gro

# Box size - copy bottom of .gro file
setmd box_l 4.09 4.09 4.09


# if gro is loaded, read it rather than creating a new system
if { $argc > 0} {
    set grofile [lindex $argv 0]
    puts "Reading $grofile"

    set fp [open $grofile r]
    set file_data [read $fp]
    close $fp

    # particle index
    set i 0
    # process data
    set data [split $file_data "\n"]
    foreach line $data {				
        if {[lindex $line 0]!=""} {
            set type 0
            # Contains: id type molecule mass pos v
            part $i pos [lindex $line 3] [lindex $line 4] [lindex $line 5] \
              v [lindex $line 6] [lindex $line 7] [lindex $line 8] \
              type $type molecule 0 mass 18.01540
            incr i
        }
    }
}


setmd time 0.0
setmd skin 0.4
setmd time_step 0.001

# Reduced unit where kT=1
thermostat langevin 1.0 1.0

# tabulated interaction 

# Note: *generating* the blockfile will require the existence of the tabulated
# file in the current directory.
inter 0 0 tabulated "table_CG_CG.tab"


# Create lists of particles for RDF
set list1 ""
for { set j 0 } { $j < [setmd n_part] } { incr j } {
    if { [part $j print type]==0 } {
        lappend list1 $j
    }
}

set out [open "| gzip -c - > conf.esp.gz" w]
blockfile $out write variable all
blockfile $out write particles [list id type molecule mass pos v]
blockfile $out write interactions
blockfile $out write thermostat
blockfile $out write tclvariable [list list1]
close $out

exit
