#!/usr/bin/env gnuplot



# Random walker over q = 1000 states



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced color size 7, 5.3 font "Helvetica, 32"
set output "walk_q1e3.eps"

reset

set style line 1 lt 1 lc rgb "#002080" lw 1.0 pt  7 ps 1.7
set style line 2 lt 2 lc rgb "#a00000" lw 1.5 pt  4 ps 1.4
set style line 3 lt 3 lc rgb "#a00040" lw 1.5 pt 12 ps 2.0
set style line 4 lt 4 lc rgb "#00a020" lw 1.5 pt  8 ps 2.0
set style line 5 lt 5 lc rgb "#a0a000" lw 1.5 pt 10 ps 2.0
set style line 7 lt 1 lc rgb "#808080"

set multiplot

# main plot
set size 1, 1
set origin 0, 0

set xtics 1 offset 0, 0.2
set mxtics 10
set xlabel "Simulation time, {/Times-Italic t} ({/Symbol \264} 10^{/*0.7 3})" offset 0, 0.5
set xrange [0:10]

set ytics 1
set mytics 10
set ylabel 'Estimated entropy, ~{/Times-Italic S}{0.5\^}' offset 1.0, 0
set yrange [4.5:7.7]

set key left bottom Left reverse spacing 1.5 at 1.6, 7.1 width 0 maxrows 2

fn = "../../data/walk/walk_q1e3.log"
ref = `head -n 1 @fn | cut -f4`
#print ref

plot [:][:] \
    fn u ($1/1e3): 2: 3 every 1 w lp ls 1 t "Uncorrected", \
    fn u ($1/1e3): 5: 6 every 1 w lp ls 2 t "Linear", \
    fn u ($1/1e3): 9:10 every 1 w lp ls 4 t "Exponential", \
    ref w l ls 7 t "Reference"


# inset
# do not reset, because we need the same line styles
set size 0.6, 0.56
set origin 0.35, 0.17

unset key

set logscale x
set logscale y

set format x "10^{/*0.7 %T}"
set xtics auto offset 0, 0.2
set mxtics default
set xlabel "Simulation time, {/Times-Italic  t}" offset 0, 0.4

set format y "10^{/*0.7 %T}"
set ytics auto
set mytics default
set ylabel '{/*1.8 |} ~{/Times-Italic S}{0.5\^} - {/Times-Italic S}^{ {/Times ref} }{/*1.8 |}' offset 0, 0

set style arrow 3 head filled size screen 0.01,15,45 ls 7
set arrow 1 from 1.4e3, 1.3 to 1e3, 0.4995 arrowstyle 3
set label 1 at 1.0e3, 2.1 "({/Times-Italic M} - 1) {/Times /} (2 {/Times-Italic t})"

plot [1e2:1e4][1e-3:6] \
    fn u ($1):(abs( $2-$4)) every 1 w lp ls 1 ps 1.0 t "Uncorrected", \
    fn u ($1):(abs( $5-$4)) every 1 w lp ls 2 ps 0.7 t "Linear", \
    fn u ($1):(abs( $9-$4)) every 1 w lp ls 4 ps 1.0 t "Exponential", \
    999/2./x ls 7 notitle

unset multiplot
unset output
set terminal pop
reset
