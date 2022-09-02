set xl 'pT'
set yl "R_{pW}"

set xr [-.2:1]
set yr [0:1.]

reaction = "{gg;qq}"
sqrs     = "18.9"
set style fill transparent solid .5

set key b
set grid 

# see figs of 1212.0434

set tit "ab -> cd -> h, pPb,  sqrt(s) = ".sqrs." TeV"
p "< paste out/RpA_\{rs=".sqrs.",A=9\}.dat out/RpA_\{rs=".sqrs.",A=184\}.dat" u 1:($2/$4) t "pT=1 GeV"

pause -1

