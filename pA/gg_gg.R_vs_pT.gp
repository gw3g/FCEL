set xl 'pT'
set yl "R_{pPb}"

set xr [0:10]
set yr [0:1.3]


set key b
set grid 

set tit "gg -> gg -> h,  sqrt(s) = 8.16 TeV"
p "out/RpA_{rs=8.16,y=0.0}.dat" u 1:3:4 w filledcurves t "y=0",\
  "out/RpA_{rs=8.16,y=2.0}.dat" u 1:2 t "y=2",\
  "out/RpA_{rs=8.16,y=4.0}.dat" u 1:2 t "y=4"

pause -1

