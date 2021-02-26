set xl '{/Symbol g}'
set yl "r"

set xr [0:4]
set yr [0.4:1.1]


set key b
set grid 

set tit "r for neutrinos"
p "out/rFCEL_gamma_{E=1.00,xi=0.5}.dat" u 1:2,\
  "out/rFCEL_gamma_{E=1.00,xi=0.5}.dat" u 1:3,\
  "out/rFCEL_gamma_{E=1.00,xi=0.5}.dat" u 1:4

pause -1

