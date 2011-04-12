set term postscript eps enhanced color
set output "microwave.eps"
set xrange [-16:16]
set yrange [-16:16]
set zrange [-2:10]
unset surface
set samples 100
set isosamples 100
set pm3d at s 
f(x)=(0.5+2.55008*x-0.588013*x**2+0.032445*x**3+0.00124411*x**4-0.0000973516*x**5)
splot f(sqrt(x**2+y**2)) title "Microwave field, 20cm x 20cm microwave"
