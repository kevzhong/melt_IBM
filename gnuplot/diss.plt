set term dumb 140 80
plot "../stringdata/diss.txt" u 1:2 w l title 'eps';
set xlabel "time"
set ylabel "epsilon"
