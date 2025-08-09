set term dumb 150 80
set logscale y
plot "../stringdata/remeshStats.txt" u 4 w l title 'dx_remesh';
set ylabel "dx_remesh"
