set term postscript eps enhanced color
set output "plot.eps"
set xrange [-0.05:0.1]
set zrange [0:600]
set cbrange [0:200]
set xlabel "x"
set ylabel "y"
set zlabel "Heat"
unset surface
set samples 100
set isosamples 100
set pm3d at s 
splot '0000407.dat' u 2:1:3  title "Bacon after 16 seconds"
