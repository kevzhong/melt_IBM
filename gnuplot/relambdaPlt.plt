set term dumb 140 80
plot "../stringdata/diss.txt" u 1:8 w l title 'relambda';
set xlabel "time"
set ylabel "relambda"
