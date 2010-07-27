# Simulation of binary mixture (explicit oil/implicit water) + extra molecule
#
#
# @author  : Tristan BEREAU (bereau@cmu.edu)
# @version : 1.00 (May 20, 2010)
# 

# Useful variables                                          #
#############################################################

set box_l    6

# Total number of oily particles
set N_oil    200
# Total number of water particles
set N_wat       0
# Total number of extra particles
set N_extra     1
set lj_sigma  1.0
set ljoffset  0.0
set wc        1.6


# Metadynamics parameters
set omega 0.001
set delta_s .5
set z_min   -6.0
set z_max    6.0
set f_bound 1.
set n_bins 200


# Integration parameters                                    #
#############################################################


setmd box_l $box_l $box_l [expr 5*$box_l]
setmd periodic 1 1 1

setmd time_step 0.001; setmd skin 0.4
set temp 1; set gamma 1
thermostat dpd $temp $gamma [expr 1.12+1.8]
#thermostat langevin $temp $gamma

# if pdb is loaded, read it rather than creating a new system
if { $argc > 0} {
    set pdbfile [lindex $argv 0]
    puts "Reading $pdbfile"

    set fp [open $pdbfile r]
    set file_data [read $fp]
    close $fp

    # particle index
    set i 1
    # process data
    set data [split $file_data "\n"]
    foreach line $data {
	if {[lindex $line 0]=="ATOM"} {
	    set type [expr int([string range [lindex $line 11] 1 end])]
	    part $i pos [lindex $line 6] \
		[lindex $line 7] \
		[lindex $line 8] \
		type $type molecule_id 0 mass 1.0
	    incr i
	}
    }
} else {
    puts "Need pdb input file."
    exit 1
}


# 0 and 1 are hydrophobic particles; 2 is inclusion
inter 0 0 lj-cos2  0.8 $lj_sigma $ljoffset $wc
inter 0 1 lj-cos2  [expr 0.7*0.8] $lj_sigma $ljoffset $wc
inter 1 1 lj-cos2  0.8 $lj_sigma $ljoffset $wc
inter 0 2 lennard-jones 1.0 $lj_sigma [expr 1.1225*$lj_sigma] [expr 0.25*1.0] 0.
inter 1 2 lennard-jones [expr 1.4*1.2] $lj_sigma [expr 1.1225*$lj_sigma] [expr 0.25*1.4*1.2] 0.

analyze set chains 1 1 $N_oil

# Define virtual particle: Center of mass of the oily phase.
part 0 pos 0 0 0 virtual 1 molecule_id 0 type 3


set part_extra_id [expr $N_oil+$N_wat+1]

metadynamics set relative_z $part_extra_id 0 $z_min $z_max $omega $delta_s $f_bound $n_bins

#  Warm up                                                  #
#############################################################
thermostat off
thermostat langevin 0.0 1.0
set warm_times     10
set warm_steps     200
puts "Warming up..."
set cap 0.01
set capincr [expr 1000./($warm_times*1.)]
for {set j 0} {$j<$warm_times} {incr j} {
    puts "Warmup run #$j out of $warm_times"
    inter ljforcecap  $cap
    integrate $warm_steps
    set cap [expr $cap+$capincr]
}
inter ljforcecap 0

for {set j 0} {$j<10} {incr j} {
    integrate 500
    puts "step $j/10"
}
setmd time_step 0.01
puts "time_step 0.01"
for {set j 0} {$j<40} {incr j} {
    integrate 500
    puts "step $j/40"
}

thermostat off
thermostat dpd $temp $gamma [expr 1.12+1.8]


setmd time 0

set out [open "|gzip -c - > conf.esp.gz" w]
blockfile $out write variable all
blockfile $out write particles [list id type molecule mass virtual pos v]
blockfile $out write interactions
blockfile $out write thermostat
close $out

exit
