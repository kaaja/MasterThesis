set terminal postscript eps color solid font "helvetica,17" linewidth 1.5
set output "forceCoeffs.eps"
set xlabel "Time"
set ylabel "Cd"
set grid
set style data line
plot \
"postProcessing/forceCoeffs/0/forceCoeffs.dat" every::400 using 1:3 title "Cd"
#    EOF
#  every ::100 (Skipping first hundred points)
