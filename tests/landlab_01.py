print("Checking Landlab...")
import landlab
grid = landlab.RasterModelGrid((4, 5), xy_spacing=(3, 4))
z = mg.add_zeros("topographic__elevation", at="node")
print("Landlab grid has %d evaluation points." % len(z))
