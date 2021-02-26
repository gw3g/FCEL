set xl 'Qs/M'
set yl "r"

set xr [0:.6]
set yr [0.4:1.1]


set key b
set grid 

set tit "r for neutrinos"
p "out/rFCEL_kappa_{E=1.00,xi=0.5}.dat" u 1:2,\
  "out/rFCEL_kappa_{E=1.00,xi=0.5}.dat" u 1:3,\
  "out/rFCEL_kappa_{E=1.00,xi=0.5}.dat" u 1:4

pause -1

