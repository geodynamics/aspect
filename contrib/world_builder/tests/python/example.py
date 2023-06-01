from gwb import WorldBuilderWrapper

import sys, getopt

def main(argv):
  filename = ""
  try:
     opts, _ = getopt.getopt(argv,"i:",["ifile="])
  except getopt.GetoptError:
     print("Failed to get commandline argument.")
     sys.exit(2)
  for opt, arg in opts:
     if opt in ("-i", "--ifile"):
        filename = arg

  world_builder = WorldBuilderWrapper(filename, False, "", 1);

  print ("2d temperature:")
  print ("temperature in Python = ", world_builder.temperature_2d(120.0e3,500.0e3,0));
  print ("3d temperature:")
  print ("temperature in Python = ", world_builder.temperature_3d(120.0e3,500.0e3,0,0));
  print ("2d composition:")
  print ("composition in Python = ", world_builder.composition_2d(120.0e3,500.0e3,0,3));
  print ("3d composition:")
  print ("composition in Python = ", world_builder.composition_3d(120.0e3,500.0e3,.0e3,0e3,3));


if __name__ == "__main__":
   main(sys.argv[1:])
