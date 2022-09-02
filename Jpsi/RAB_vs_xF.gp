set xl 'x_F'
set yl "R_{A/B}(x_F)"

set xr [-.2:1]
set yr [0:1.]

sqrs     = "38.7"

set key b
set key l
set grid 

# see figs of 1212.0434

set tit "J/psi,  sqrt(s) = ".sqrs." GeV"
p "out/R_WBe_{rs=38.7,pt=1.0}.dat" u 1:4 t "W/Be, pT=1 GeV",\
  "out/R_FeBe_{rs=38.7,pt=1.0}.dat" u 1:4 t "Fe/Be, pT=1 GeV"

pause -1

