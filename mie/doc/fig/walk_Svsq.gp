#!/usr/bin/env gnuplot



# Potts model



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced color size 7, 5 font "Helvetica, 36"
set output "walk_Svsq.eps"

reset

set style line 1 lt 1 lc rgb "#002080" lw 1.0 pt  7 ps 1.7
set style line 2 lt 2 lc rgb "#a00000" lw 1.5 pt  4 ps 1.4
set style line 3 lt 3 lc rgb "#a00040" lw 1.5 pt 12 ps 2.0
set style line 4 lt 4 lc rgb "#00a020" lw 1.5 pt  8 ps 2.0
set style line 5 lt 5 lc rgb "#a0a000" lw 1.5 pt 10 ps 2.0
set style line 7 lt 1 lc rgb "#404040"

set logscale x
set format x "10^{/*0.7 %T}"
#set xtics 1 offset 0, 0.2
#set mxtics 2
set xlabel "Number of states, {/Times-Italic M}" offset 0, 0.
#set xrange [0:10]

set ytics 1
set mytics 2
set ylabel 'Estimated entropy, ~{/Times-Italic S}{0.5\^}' offset 1.0, 0
set yrange [4:14]

set key left bottom Left reverse spacing 1.5 at 1e2, 11.8 width -2 maxrows 2

fn = "../../data/walk/walk_Svsq.dat"
#print ref

plot [:][:] \
    fn u ($1): 2: 3 every 1 w lp ls 1 ps 3   t "Uncorrected", \
    fn u ($1): 5: 6 every 1 w lp ls 2 ps 2   t "Linear", \
    fn u ($1): 9:10 every 1 w lp ls 4 ps 3   t "Exponential", \
    fn u ($1):4  w l ls 7 t "Reference"


unset output
set terminal pop
reset
