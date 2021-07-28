set xl 'E/GeV'
set xr [1e3:1e10]
set yr [0.7:0.9]
set log x


set key t
set grid 
reaction = "{gg;qq}"
kTval    = "2.00"
x2val    = "1e-05"

set tit "r for neutrinos: k_T=".kTval.", x_2=".x2val." (in legend: d{/Symbol s}/dx_F; CR flux)"
p "out/r2_nu_FCEL_GSF_".reaction."_{kT=".kTval.",x2=1e-05}.dat"  u 1:4 t " AJP; Dembinski ",\
  "out/r2_nu_FCEL_H3a_".reaction."_{kT=".kTval.",x2=1e-05}.dat"  u 1:4 t " AJP; H3a ",\
  "out/r2_nu_FCEL_knee_".reaction."_{kT=".kTval.",x2=1e-05}.dat" u 1:4 t " AJP; broken-power " ,\
  "out/r2_nu_scaling1_".reaction."_{kT=".kTval.",x2=1e-05}.dat"  u 1:4 w lines t " scaling - n=6; {/Symbol f} = E^{-2.7} ",\
  "out/r2_nu_scaling2_".reaction."_{kT=".kTval.",x2=1e-05}.dat"  u 1:4 w lines t " scaling - n=6; {/Symbol f} = E^{-3} ",\
  "out/r2_nu_scaling3_".reaction."_{kT=".kTval.",x2=1e-05}.dat"  u 1:4 w lines t " scaling - n=6; {/Symbol f} = E^{-3.4} "

pause -1

