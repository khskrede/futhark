#! /bin/bash
#This file animates all the png files
file_out=animate.avi
fps=24
obr=1200
opt="mbd=2:keyint=132:cmp=2:subcmp=2:dia=2:mv0:last_pred=3"
codec="msmpeg4v2"

mencoder -ovc lavc -lavcopts vcodec=$codec:vbitrate=$obr:vpass=1:$opt:turbo \
-nosound -mf fps=$fps:type=png   -o NUL: mf://data/*.dat.png
mencoder -ovc lavc -lavcopts vcodec=$codec:vbitrate=$obr:vpass=2:$opt \
-nosound -mf fps=$fps:type=png -o $file_out mf://data/*.dat.png
