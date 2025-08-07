set term dumb 140 80
plot "../stringdata/diss.txt" u 1:3 w l title 'tke';
set xlabel "time"
set ylabel "tke"
