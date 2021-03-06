set xl 'pT'
set yl "R_{pPb}"

set xr [-6:6]
set yr [0:1.3]


set key b
set grid 

# see fig.2 of 2003.06337

set tit "gg -> gg -> h,  sqrt(s) = 8.16 TeV"
p "out/RpA_{rs=8.16,pT=2.0}.dat" u 1:3:4 w filledcurves t "pT=2 GeV",\
  "out/RpA_{rs=8.16,pT=6.0}.dat" u 1:3:4 w filledcurves t "pT=6 GeV"

pause -1

