#!/bin/bash

for file in ./data/*.dat
do
	echo "set title \"3D data (2D Heat Map)\"
	set xlabel \"x\"
	set ylabel \"y\"
	set zlabel \"heat\"
	set terminal png size 800, 600 enhanced
	set output \"${file}.png\"
	set cbrange [0:400]
	set zrange [0:400]
	splot \"${file}\" using 2:1:3 with pm3d " | gnuplot
done

