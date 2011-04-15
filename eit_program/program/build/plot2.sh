#!/bin/bash

for file in ./data/*.flo
do
	echo "set title \"3D data (2D Heat Map)\"
	set xlabel \"x\"
	set ylabel \"y\"
	set zlabel \"flow\"
	set terminal png size 800, 600 enhanced
	set output \"${file}.png\"
	#set cbrange [10.0:10.0]
	#set zrange [10.0:10.0]
	splot \"${file}\" using 2:1:3 with pm3d " | gnuplot
done

