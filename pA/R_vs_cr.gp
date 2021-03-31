set xl 'C_R'
set yl "R_{pPb}"

set xr [0:20]
set yr [0:1.]

reaction = "{gg;gg}"
sqrs     = "8.16"

set key b
set grid 

set label "8" at 3-.1, graph .32
set arrow from 3, graph .3 to 3, graph 0

set label "3" at 4./3-.1, graph .32
set arrow from 4./3, graph .3 to 4./3, graph 0

set label "6" at 10./3-.1, graph .32
set arrow from 10./3, graph .3 to 10./3, graph 0

set label "10" at 6-.2, graph .32
set arrow from 6, graph .3 to 6, graph 0

set label "15" at 16./3-.2, graph .32
set arrow from 16./3, graph .3 to 16./3, graph 0

set label "27," at 8.-.3, graph .32
set arrow from 8, graph .3 to 8, graph 0

set label "24" at 25./3, graph .32
set arrow from 25./3, graph .3 to 25./3, graph 0

set label "15" at 28./3-.2, graph .32
set arrow from 28./3, graph .3 to 28./3, graph 0

set label "21" at 40./3-.2, graph .32
set arrow from 40./3, graph .3 to 40./3, graph 0

set label "35" at 12-.2, graph .32
set arrow from 12, graph .3 to 12, graph 0

set label "42" at 34./3-.2, graph .32
set arrow from 34./3, graph .3 to 34./3, graph 0

# see figs of 2003.06337

set tit "ab -> cd -> h, pPb,  sqrt(s) = ".sqrs." TeV"
p "out/R_CR_{rs=".sqrs.",pT=0.0,y=-3.5}.dat" u 1:3:4 w filledcurves t "pT=0 GeV, y=-3.5",\
  "out/R_CR_{rs=".sqrs.",pT=0.0,y=3.5}.dat" u 1:3:4 w filledcurves t "pT=0 GeV, y=+3.5",\
  "out/R_CR_{rs=".sqrs.",pT=2.0,y=-3.5}.dat" u 1:3:4 w filledcurves t "pT=2 GeV, y=-3.5",\
  "out/R_CR_{rs=".sqrs.",pT=2.0,y=3.5}.dat" u 1:3:4 w filledcurves t "pT=2 GeV, y=+3.5"

pause -1

