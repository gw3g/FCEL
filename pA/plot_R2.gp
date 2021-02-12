set xl 'y'
set yl "R_{pPb}"

set xr [-6:6]
set yr [0:1.3]


set key b
set grid 



set tit "gg -> gg -> h,  sqrt(s) = 8.16 TeV"
p "R_pT2GeV.dat" u 1:3:4 w filledcurves t "pT=2",\
  "R_pT6GeV.dat" u 1:3:4 w filledcurves t "pT=6"

pause -1

