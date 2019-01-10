#!/bin/sh

../bin/tabulate 3 0.501 0.75 0.001 > tabulate_3.dat
../bin/tabulate 4 0.501 0.75 0.001 > tabulate_4.dat
gnuplot tabulate.gp
gvfs-open tabulate.png
