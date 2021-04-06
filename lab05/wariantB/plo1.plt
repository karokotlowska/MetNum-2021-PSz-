set term jpeg size 800,600

set xlabel 'm-ilosc iteracji'
set ylabel 'lambda'

set out 'wynik_lin.png'

#set linetype 1 lc rgb "red"

set xtics 1
set ytics 2

set yrange [0:4]

plot 'lambda.dat' u 1:2 w l t 'przyblizenia wartosci wlasnych'