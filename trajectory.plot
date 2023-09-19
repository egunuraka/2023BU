do for [ii=1:365]{
	splot "trajectory" u 2:3:4 every ::::ii w linespoint t "earth", \
	"trajectory" u 5:6:7  every ::::ii w linespoint t "2023BU", \
	"trajectory" u 8:9:10 every ::::ii w points t "sun", \
	"trajectory" u 11:12:13 every ::::ii w linespoint t "moon"
	pause 0.5 
}
pause -1
