set term dumb 150 80
stats '../stringdata/phasefield.txt' using 2 every ::0::0 nooutput
v0 = STATS_min
plot "../stringdata/phasefield.txt" u 1:($2/v0) w l title 'vol';
set xlabel "time"
set ylabel "vol"
