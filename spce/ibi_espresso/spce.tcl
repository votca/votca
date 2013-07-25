# check if required features are compiled in:

proc require_feature {feature} {
    global errf
    if { ! [regexp $feature [code_info]]} {
	puts "not compiled in: $feature"
	exit -42
    }
}
#write output file 
proc write_data {file_name i time_step} {
    set f [open $file_name "a"]
    puts $f " [expr $i*$time_step]  [expr [analyze energy total]/[setmd n_part]] [expr (([analyze energy total] - [analyze energy kinetic])/[setmd n_part])] [expr [analyze energy kinetic]/[setmd n_part]] "
    close $f
}

proc write_xyz {file_name} {
  set f [open $file_name "a"]
  puts $f "[setmd n_part]"
  puts $f "xyz by Espresso, time [setmd time]"
  for { set i 1} { $i <= [setmd n_part] } { incr i 1} {
    puts $f "CG [part $i print pos]"
  }
  close $f
}

#check for necessary feature
require_feature "TABULATED"
require_feature "MASS"
#write code info to terminal
puts "Program Information: \n[code_info]\n"


# set system properties:

set box_length 4.053
set skin 0.5
set temperature 2.49 
set time_step 0.001
set friction 5.0
set int_steps 1000
set int_steps 1000

setmd box_l $box_length $box_length $box_length
setmd skin $skin
setmd time_step $time_step
thermostat langevin $temperature $friction

#open file to read
set infile [open "spce.esp" r]

#skip comment lines
##comment out first two lines in .esp, when converted from .gro
set in_string "#"
while { [string index $in_string 0] == "#" } { gets $infile in_string }

#read data lines
while { 1 } {
  # if at file end, close file and stop loop
  if {[eof $infile]} {
    close $infile
    break
  }
  # write numbers into variables
  scan $in_string "%s %s %i %f %f %f %f %f %f" p_name p_sp p_id pos_x pos_y pos_z v_x v_y v_z
  # determine type (all particles of same type)
  set p_type 0
  # set up particle
  part $p_id pos $pos_x $pos_y $pos_z type $p_type v $v_x $v_y $v_z mass 18.01540
  # read next line
  gets $infile in_string
}

inter $p_type $p_type tabulated "CG_CG.tab" 

integrate 0
#measure simulation time
set starttime [ clock clicks -milliseconds ]
#warmup loop
puts "number of particles [setmd n_part]"
set temperature [expr 2.0/3.0*[analyze energy kinetic]/[setmd n_part]]
puts "inital temperature: [format "%.2f" $temperature]"
for { set i 0} { $i < 100 } { incr i 1} {
  puts -nonewline "warmup step $i \r"
  flush stdout;
  integrate 10
}
set temperature [expr 2.0/3.0*[analyze energy kinetic]/[setmd n_part]]
puts "Running at temperature T=[format "%.2f" $temperature]"

set f [open "traj.xyz" "w"]
close $f
set f [open "energy.dat" "w"]
close $f

#integration loop
for { set i 0} { $i < $int_steps } { incr i 1} {
  #puts "integration step $i "
  #10 int steps on "C-level" for 1 int step on "TCL level"
  integrate 10
  flush stdout;
  #write some output to terminal
  puts -nonewline "time: [format "%.3f" [ expr $i*$time_step ]] potential energy: [format "%.2f" [ expr (([analyze energy total] - [analyze energy kinetic]) / [setmd n_part]) ]]\r "
  #store energies in energy.dat
  write_data "energy.dat" $i $time_step
  #for paraview: usage open files -> part_data* -> use glyph ->spheres radius 0.05
  #writevtk part_data$i.vtk
  #for vmd output
  write_xyz traj.xyz
  #imd positions
}
puts "\n"

set stoptime [ clock clicks -milliseconds ]

puts "time per MD step: [ expr (($stoptime-$starttime)/$int_steps)] milliseconds"
puts "pot. energy: [format "%.2f" [ expr (([analyze energy total] - [analyze energy kinetic]) / [setmd n_part]) ]]"
puts "final temperature: [format "%.2f" [expr 2.0/3.0*[analyze energy kinetic]/[setmd n_part] ] ]"
puts "Finished."





