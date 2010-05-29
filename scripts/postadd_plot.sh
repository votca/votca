do_external postadd dummy "$1" "$2"

killall gnuplot_x11

gnuplot -geometry 1024x700 -persist ../allplots.gp 

true
