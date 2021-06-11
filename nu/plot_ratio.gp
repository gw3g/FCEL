set xl 'E/GeV'
set xr [1e3:1e10]
set yr [0.7:0.9]
set log x


set key t
set grid 

set tit "r for neutrinos: k_T=0, x_2=10^{-5} (in legend: d{/Symbol s}/dx_F; CR flux)"
p "out/r_nu_FCEL_GSF_{qg;qg}_{kT=0.00,x2=1e-05}.dat"  u 1:4 t " AJP; Dembinski ",\
  "out/r_nu_FCEL_H3a_{qg;qg}_{kT=0.00,x2=1e-05}.dat"  u 1:4 t " AJP; H3a ",\
  "out/r_nu_FCEL_knee_{qg;qg}_{kT=0.00,x2=1e-05}.dat" u 1:4 t " AJP; broken-power " ,\
  "out/r_nu_scaling1_{qg;qg}_{kT=0.00,x2=1e-05}.dat"  u 1:4 w lines t " scaling - n=6; {/Symbol f} = E^{-2.7} ",\
  "out/r_nu_scaling2_{qg;qg}_{kT=0.00,x2=1e-05}.dat"  u 1:4 w lines t " scaling - n=6; {/Symbol f} = E^{-3} ",\
  "out/r_nu_scaling3_{qg;qg}_{kT=0.00,x2=1e-05}.dat"  u 1:4 w lines t " scaling - n=6; {/Symbol f} = E^{-3.4} "

pause -1

