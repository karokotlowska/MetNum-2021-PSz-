set term png # ustawienie typu terminala, np. m.in. x11 (ekran), postscript, pdf, png, table (kolumny współrzędnych). 

set out "rys.png" # ustawienie nazwy pliku wyjściowego 
set ytics nomirror        # Uniezaleznienie znacznikow skali osi OY od drugiej osi OY
set ylabel "||r_k||_2" tc rgb "red" offset 2, 0      # Podpis LEWEJ osi OY 
set y2tics nomirror       # Analogicznie.
set logscale y            # Skala logarytmiczna na LEWEJ osi OY
set xl "k - numer iteracji" # tytuł osi x
set yl "||r_k||_2" # tytuł osi y
set title "Wychylenie x(t)" # tytuł wykresu
set ytics nomirror 

plot "result.dat" u 1:2 w lp ls 1 t "||r_k||_2,", \
"" u 1:5 w lp ls 2 t "||x_k||_2" axes x1y2

#plot "result.dat" u 1:2 w lp ls 1 t "||r_k||_2,", \
"" u 1:5 w lp ls 2 t "||x_k||_2" axes x1y2

#p "result.dat" u 1:2 w p lp 3 pt 6 t "h=0.1" # rysowanie wykresu
