#!/bin/sh

# consider using ffmpeg, it's simpler than this.

#mencoder "mf://$1/*.png" -mf fps=$3 -o $2 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=9600
#mencoder "mf://$1/*.png" -mf fps=$3 -o $2 -ovc lavc -lavcopts vcodec=msmpeg4v2
#mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc raw -oac copy -o output.avi
#mencoder mf://*.png -mf w=1275:h=1275:fps=$3:type=png -ovc raw -oac copy -o $2
#WORKS FOR ALMOST EVERYTHING
#mencoder mf://*.png -mf w=1275:h=1275:fps=$3:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $2
#WHEN THERE ARE QUALITY ISSUES
#mencoder mf://*.png -mf w=1275:h=1275:fps=$3:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:vbitrate=5000 -oac copy -o $2
mencoder mf://$1/*.png -mf w=1275:h=1275:fps=$3:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:vbitrate=5000 -oac copy -o $2
