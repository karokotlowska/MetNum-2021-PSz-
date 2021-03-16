set terminal png

set output "f_r.png"
set autoscale
set logscale y
set xlabel "k"
set ylabel "||r_k||_2"
set title "||r_k||_2 = f(k)"  
plot "a.dat" u 1:2 w l t "f(k)",\
	 "b.dat" u 1:2 w l t "f(k)"


set output "f_x.png"

set xlabel "k"
set ylabel "||x_k||_2"
set yrange [3635:3665]
set title "||x_k||_2 = f(k)"  
plot "a.dat" u 1:4 w l t "f(k)",\
    "b.dat" u 1:4 w l t "f(k)"




