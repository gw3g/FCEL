set xl 'pT'
set yl "R_{pPb}"

set xr [-6:6]
set yr [0:1.3]

reaction = "{qg;qg}"
sqrs     = "8.16"
set style fill transparent solid .5

set key b
set grid 

# see figs of 2003.06337

set tit "ab -> cd -> h, pPb,  sqrt(s) = ".sqrs." TeV"
p "out/RpA_".reaction."_{rs=".sqrs.",pT=2.0}.dat" u 1:3:4 w filledcurves t "pT=0 GeV",\
  "out/RpA_".reaction."_{rs=".sqrs.",pT=6.0}.dat" u 1:3:4 w filledcurves t "pT=2 GeV"

pause -1

