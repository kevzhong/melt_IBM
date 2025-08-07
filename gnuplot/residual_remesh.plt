set term dumb 150 80
set logscale y
plot "../stringdata/remeshStats.txt" u 0:(abs($3)) w l title 'residual';
set ylabel "residual"