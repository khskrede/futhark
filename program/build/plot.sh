#!/bin/bash

for (( i=0; i<$1; i++ )) 
do
	echo "set title \"3D data (2D Heat Map)\"
	set xlabel \"x\"
	set ylabel \"y\"
	set zlabel \"heat\"
	set terminal png size 800, 600 enhanced
	set output \"data/${i}.png\"
	set cbrange [0:200]
	set zrange [0:200]
	splot \"data/${i}.dat\" using 1:2:3 with pm3d " | gnuplot
done

rm ./data/*.dat
ffmpeg -f image2 -i ./data/%d.png ./data/heat.mpg
rm ./data/*.png
