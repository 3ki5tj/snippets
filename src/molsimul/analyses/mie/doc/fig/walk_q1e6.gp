#!/usr/bin/env gnuplot



# Potts model



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced color size 9, 7 font "Helvetica, 36"
set output "walk_q1e6.eps"

reset

set style line 1 lt 1 lc rgb "#002080" lw 1.0 pt  7 ps 1.7
set style line 2 lt 2 lc rgb "#a00000" lw 1.5 pt  4 ps 1.4
set style line 3 lt 3 lc rgb "#a00040" lw 1.5 pt 12 ps 2.0
set style line 4 lt 4 lc rgb "#00a020" lw 1.5 pt  8 ps 2.0
set style line 5 lt 5 lc rgb "#a0a000" lw 1.5 pt 10 ps 2.0
set style line 7 lt 1 lc rgb "#808080"

set multiplot

set size 1, 1
set origin 0, 0

set xtics 1 offset 0, 0.2
set mxtics 2
set xlabel "Simulation time, {/Times-Italic t} ({/Symbol \264} 10^{/*0.7 3})" offset 0, 0.5
set xrange [0:10]

set ytics 1
set mytics 2
set ylabel 'Estimated entropy, ~{/Times-Italic S}{0.5\^}' offset 1.0, 0
set yrange [4:14]

set key left bottom Left reverse spacing 1.5 at 1.0, 4.5 width -8 maxrows 3

fn = "../../data/walk/walk_q1e6.log"
ref = `head -n 1 @fn | cut -f4`
#print ref

plot [:][:] \
    fn u ($1/1e3): 2: 3 every 1 w lp ls 1 pt  7 ps 2   t "Uncorrected", \
    fn u ($1/1e3): 5: 6 every 1 w lp ls 2 pt  4 ps 1.4 t "Linear, 1st order", \
    fn u ($1/1e3): 7: 8 every 1 w lp ls 3 pt 12 ps 2   t "Linear, 2nd order", \
    fn u ($1/1e3): 9:10 every 1 w lp ls 4 pt  8 ps 2   t "Exponential,  1st order", \
    fn u ($1/1e3):11:12 every 1 w lp ls 5 pt 10 ps 2   t "Exponential,  2nd order", \
    ref w l ls 7 t "Reference"


unset multiplot
unset output
set terminal pop
reset
