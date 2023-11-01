## This is a script written by Thilo Wrona which rewrites FastScape VTK files
## into VTS files that include the correct time of ASPECT solution files. This is
## necessary because the VTK files do not include a time and cannot be viewed
## simultaneously with ASPECT in paraview.

import sys
import os
import glob
import numpy as np

import pandas as pd
from xml.etree import ElementTree as ET

#### import the simple module from the paraview, change path.
sys.path.insert(1,'/path/to/paraview/lib/python3.8/site-packages/')
from paraview.simple import *


def get_times_pvd(filename):
    # Read in the solution.pvd file
    data = pd.read_csv(filename,header=5,delim_whitespace=True, names=[0, 1, 2, 3, 4])
    data = data.to_numpy()

     # Remove text characters and the final two rows which don't contain any time info.
    times = np.zeros(len(data[:,1])-2)
    for i in range(data.shape[0]-2):
        temp_str = data[i,1]

        # Column 1 contains a string of " "timestep=0" "
        # Here we remove the first ten characters, and the final character.
        temp_str = temp_str[10:]
        temp_str = temp_str[:-1]

        times[i] = float(temp_str)

    return times 


#%% Get file paths (absolute, not relative paths!)
def absoluteFilePaths(directory):
   for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           yield os.path.abspath(os.path.join(dirpath, f))

# This will work if you run the script inside the ASPECT output folder,
# otherwise you can put the absolute path.
# This will make sure we only pick up the topography files.
FileNames = glob.glob('./VTK/Topography*.vtk')
FileNames = sorted(FileNames)

DirPath = os.path.dirname(os.path.realpath(__file__))

#%%

# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

for File in FileNames:

    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()
    
    # create a new 'Legacy VTK Reader'
    topography0000 = LegacyVTKReader(FileNames=File)
    
    # get animation scene
    animationScene1 = GetAnimationScene()
    
    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()
    
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [882, 554]
    
    # get layout
    layout1 = GetLayout()
    
    # show data in view
    topography0000Display = Show(topography0000, renderView1, 'StructuredGridRepresentation')
    
    # trace defaults for the display properties.
    topography0000Display.Representation = 'Surface'
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    # get the material library
    materialLibrary1 = GetMaterialLibrary()
    
    # show color bar/color legend
    topography0000Display.SetScalarBarVisibility(renderView1, True)
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # get color transfer function/color map for 'H'
    hLUT = GetColorTransferFunction('H')
    
    # get opacity transfer function/opacity map for 'H'
    hPWF = GetOpacityTransferFunction('H')
    
    # save data
    SaveData(File[:-4] + '.vts', proxy=topography0000, ChooseArraysToWrite=1,
        PointDataArrays=['HHHHH','basement','catchment','drainage_area','erosion_rate','topography','total_erosion'])
    
    #### saving camera placements for all active views
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [225625.0, 25625.0, 878497.0779980461]
    renderView1.CameraFocalPoint = [225625.0, 25625.0, 1135.7277145385742]
    renderView1.CameraParallelScale = 227077.82689023562
    
    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    
    Delete(topography0000)
    del topography0000

#%% Get time steps
times = get_times_pvd('solution.pvd')

#%% Write pvd file

root_attributes = {"type": "Collection", "version": "0.1", "ByteOrder": "LittleEndian"}

root = ET.Element("VTKFile", attrib=root_attributes)
root.text="\n"


doc = ET.SubElement(root, "Collection")
doc.text="\n"

for n, File in enumerate(FileNames):
    dataset_attributes = {"timestep": str(times[n]), "group": "", "group": "", "file": os.path.basename(File)[:-4]+'.vts'}
    
    ET.SubElement(doc, "DataSet", attrib=dataset_attributes).text="\n"


tree = ET.ElementTree(root)
tree.write('./VTK/topography.pvd', xml_declaration=True, encoding='utf-8', method="xml")   
    
    
    





