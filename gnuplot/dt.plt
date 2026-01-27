set term dumb 140 80
plot "../stringdata/dt.txt" u 1:2 w l title 'dt';
set xlabel "time"
set ylabel "dt"
