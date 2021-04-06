#niewykorzystywany tu
set terminal post enhanced colour solid font 20  # wybor formatu, w jakim chcemy utworzyc wynikowy rysunek

#set terminal post enhanced colour

set output "wynik.png" # nazwa pliku wynikowego

#wartosc wlasna

set title "Wartosci wlasne" # tytul wykresu
set xlabel "Numer iteracji k" # etykieta osi OX
set ylabel "{/Symbol l}_k" # etykieta osi OY

set grid
plot [1 to 300] [0 to 1.5] "przyblizeniaLambda.dat" u 1:3 w l lc 2 t "{/Symbol l}_2", "przyblizeniaLambda.dat" u 1:4 w l lc 3 t "{/Symbol l}_3", "przyblizeniaLambda.dat" u 1:5 w l lc 4 t "{/Symbol l}_4", "przyblizeniaLambda.dat" u 1:6 w l lc 5 t "{/Symbol l}_5", "przyblizeniaLambda.dat" u 1:7 w l lc 6 t "{/Symbol l}_6"

plot [1 to 300] [3.588 to 3.598] "przyblizeniaLambda.dat" u 1:2 w l lc 1 t "{/Symbol l}_1"