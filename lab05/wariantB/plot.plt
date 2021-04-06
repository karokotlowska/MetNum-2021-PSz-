#set datafile separator ","
set term png
set output "plotw.png"
set yrange [0:4]

plot 'lambda.dat' u 1:2 w p lt 3 lw 2 t 'Przybliżenia wartości własnych'