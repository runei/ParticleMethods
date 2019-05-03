# !/bin/bash

cd output
arq=`ls -v`

for a in $arq
do
	echo "set xrange[-1:1]" > input.gnu
	echo "set yrange[-1:1]" >> input.gnu
	echo "set key off" >> input.gnu
	# echo "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb 'black' behind " >> input.gnu
	# echo "set palette" >> input.gnu
	# echo "plot '"$a"' using 1:2 with points pt 7 ps variable" >> input.gnu
	echo "plot '"$a"' using 1:2:3 with points pt 7 ps variable" >> input.gnu
	echo "pause 1" >> input.gnu
	gnuplot < input.gnu
done
rm input.gnu
