set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
set output "compare.png"
p'centerline_OFC.dat' u 3:($6/1000) w lp ,'fort.11' u 3:($5/1000) w lp,'fort.11' u 3:($5/1000) w lp
