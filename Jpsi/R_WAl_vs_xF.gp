set xl 'x_F'
set yl "R_{A/B}(x_F)"

set xr [-.2:1]
set yr [0:1.]

sqrs     = "18.9"

set key b
set key l
set grid 

set tit "J/psi,  sqrt(s) = ".sqrs." GeV"
p "out/R_WAl_{rs=18.9,pt=0.0}.dat" u 1:4 t "W/Be, pT=0 GeV",\
  "out/R_WAl_{rs=18.9,pt=1.0}.dat" u 1:4 t "W/Be, pT=1 GeV",\
  "out/R_WAl_{rs=18.9,pt=2.0}.dat" u 1:4 t "W/Be, pT=2 GeV",\
  "out/R_WAl_{rs=18.9,pt=3.0}.dat" u 1:4 t "W/Be, pT=3 GeV",\
  "out/R_WAl_{rs=18.9,pt=4.0}.dat" u 1:4 t "W/Be, pT=4 GeV"

pause -1

