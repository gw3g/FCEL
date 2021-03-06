set xl 'pT'
set yl "R_{pPb}"

set xr [0:10]
set yr [0:1.3]


set key b
set grid 

# see fig.3 of 2003.06337

set tit "gg -> gg -> h,  sqrt(s) = 8.16 TeV"
p "out/RpA_{rs=8.16,y=0.0}.dat" u 1:3:4 w filledcurves t "y=0",\
  "out/RpA_{rs=8.16,y=3.0}.dat" u 1:2 t "y=3",\
  "out/RpA_{rs=8.16,y=5.0}.dat" u 1:2 t "y=5"

pause -1

