set term dumb 140 80
plot "../stringdata/part_vel.txt" u 1:3 w l title 'vx', \
	"../stringdata/part_vel.txt" u 1:4 w l title 'vy',\
	"../stringdata/part_vel.txt" u 1:5 w l title 'vz';
