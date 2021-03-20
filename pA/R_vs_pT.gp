set xl 'pT'
set yl "R_{pPb}"

set xr [0:10]
set yr [0:1.3]

reaction = "{gg;gg}"
sqrs     = "8.16"

set style fill transparent solid .5

set key b
set grid 

# see figs of 2003.06337

set tit "ab -> cd -> h,  pPb, sqrt(s) = ".sqrs." TeV"
p "out/RpA_".reaction."_{rs=".sqrs.",y=0.0}.dat" u 1:3:4 w filledcurves t "y=0",\
  "out/RpA_".reaction."_{rs=".sqrs.",y=3.0}.dat" u 1:3:4 w filledcurves t "y=3",\
  "out/RpA_".reaction."_{rs=".sqrs.",y=5.0}.dat" u 1:3:4 w filledcurves t "y=5"

pause -1

