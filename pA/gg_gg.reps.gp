set xl 'y'
set yl "R_{pPb}"

set xr [-6:6]
set yr [0:1.3]


set key b
set grid 



set tit "gg -> gg -> h,  sqrt(s) = 8.16 TeV"
p "out/R_reps_{rs=8.16,pT=2.0}.dat" u 1:2 t "1",\
  "out/R_reps_{rs=8.16,pT=2.0}.dat" u 1:3 t "8",\
  "out/R_reps_{rs=8.16,pT=2.0}.dat" u 1:4 t "27",\
  "out/R_reps_{rs=8.16,pT=2.0}.dat" u 1:5 t "ave"

pause -1

