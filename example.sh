
./state.out 10 5 5 1.0 1.0 4.0 2.0 0.0 0.1 0.5 1


cat > plot.gnu << EOF
set term  post colo dashed eps enhanced "Helvetica, 20"
set out "fig.eps"
set size 1.0,1.0
set ylabel "J"
set xrange[0:50]
set xlabel "t"

plot 'fin_rand_monitor-10-5-5-1.00-1.00-4.00-2.00-0.0000-0.1000.long' u 1:2 w l lw 2.5 lt 1 lc 0 t 'U=4.0,V=2.0, long',\
     'fin_rand_monitor-10-5-5-1.00-1.00-4.00-2.00-0.0000-0.1000.out' u 1:2 w l lw 3.5 lt 1 lc 1 t 'U=4.0,V=2.0,short'

EOF

gnuplot plot.gnu

gv fig.eps
