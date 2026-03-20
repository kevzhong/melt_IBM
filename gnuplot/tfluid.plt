set term dumb 150 80
plot "../stringdata/phasefield.txt" u 1:3 w l title 'tfluid';
set xlabel "time"
set ylabel "tfluid"
