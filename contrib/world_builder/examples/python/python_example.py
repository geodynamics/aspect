from gwb import WorldBuilderWrapper

filename = "../../tests/data/continental_plate.wb"
world_builder = WorldBuilderWrapper(filename);

print ("2d temperature:")
print ("temperature in Python = ", world_builder.temperature_2d(120.0e3,500.0e3,0));
print ("3d temperature:")
print ("temperature in Python = ", world_builder.temperature_3d(120.0e3,500.0e3,0,0));
print ("2d composition:")
print ("composition in Python = ", world_builder.composition_2d(120.0e3,500.0e3,0,3));
print ("3d composition:")
print ("composition in Python = ", world_builder.composition_3d(120.0e3,500.0e3,.0e3,0e3,3));
