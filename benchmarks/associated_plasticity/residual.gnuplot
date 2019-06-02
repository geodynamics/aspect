
!grep "Relative nonlinear" output_implicit/log.txt >output_implicit/residual
!grep "Relative nonlinear" output_explicit/log.txt >output_explicit/residual
!grep "Relative nonlinear" output_disabled/log.txt >output_disabled/residual

set log y

plot "output_implicit/residual" using 10 w lp title "implicit", \
"output_explicit/residual" using 10 w lp title "explicit", \
"output_disabled/residual" using 10 w lp title "disabled"

pause -1
