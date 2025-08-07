set term dumb 150 80
plot "../stringdata/remeshStats.txt" u 0:5 w l title 'Nrefresh';
set ylabel "Nrefresh"