set term dumb 150 80
plot "../stringdata/part_vel.txt" u 1:2 w l title 'dt';
set xlabel "time"
set ylabel "dt"
