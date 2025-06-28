set term dumb 150 80
stats '../stringdata/part_vol.txt' using 2 every ::0::0 nooutput
v0 = STATS_min
plot "../stringdata/part_vol.txt" u 1:($2/v0) w l title 'vol';
set xlabel "time"
set ylabel "vol"
