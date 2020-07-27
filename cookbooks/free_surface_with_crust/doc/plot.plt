set term png font arial 16
set output "topography.png"

plot "../output-free_surface_with_crust/statistics" using 2:15 with lines lw 2 title "Max topo crust [m]", \
     "../output-free_surface_with_crust/statistics" using 2:14 with lines lw 2 title "Min topo crust [m]", \
     "../../free_surface/output-free_surface/statistics" using 2:15 with lines lw 2 title "Max topo no crust [m]", \
     "../../free_surface/output-free_surface/statistics" using 2:14 with lines lw 2 title "Min topo no crust [m]"
