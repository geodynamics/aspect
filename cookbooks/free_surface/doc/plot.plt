set term png font arial 16
set output "free_surface_topography.png"

plot "../output-free_surface/statistics" using 2:15 with lines lw 2 title "Maximum topography [m]"
