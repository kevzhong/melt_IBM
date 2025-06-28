set term dumb 140 80
plot "../stringdata/part_pos.txt" u 1:3 w l title 'x', \
	"../stringdata/part_pos.txt" u 1:4 w l title 'y',\
	"../stringdata/part_pos.txt" u 1:5 w l title 'z';
