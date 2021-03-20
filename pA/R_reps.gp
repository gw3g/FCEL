set xl 'y'
set yl "R_{pPb}"

set xr [-6:6]
set yr [0:1.3]

reaction = "{gg;qq}"
sqrs     = "8.16"

set key b
set grid 

# see figs of 2003.06337

set tit "ab -> cd -> h, pPb, sqrt(s) = ".sqrs." TeV"
p "out/R_reps_".reaction."_{rs=".sqrs.",pT=2.0}.dat" u 1:2 t "1",\
  "out/R_reps_".reaction."_{rs=".sqrs.",pT=2.0}.dat" u 1:3 t "8",\
  "out/R_reps_".reaction."_{rs=".sqrs.",pT=2.0}.dat" u 1:4 t "27",\
  "out/R_reps_".reaction."_{rs=".sqrs.",pT=2.0}.dat" u 1:5 t "ave"

pause -1

