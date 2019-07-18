#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Sphere'
sphere1 = Sphere()

# Properties modified on sphere1
sphere1.Center = [0.5, 0.5, 0.5]
sphere1.Radius = 0.15
sphere1.ThetaResolution = 64 
sphere1.PhiResolution = 64 

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1118, 799]

# show data in view
sphere1Display = Show(sphere1, renderView1)
# trace defaults for the display properties.
sphere1Display.AmbientColor = [0.0, 0.0, 0.0]
sphere1Display.ColorArrayName = [None, '']
sphere1Display.EdgeColor = [0.0, 0.0, 0.0]
sphere1Display.CubeAxesColor = [0.0, 0.0, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# change representation type
sphere1Display.SetRepresentationType('Surface With Edges')

# save data
SaveData('./sphere.stl', proxy=sphere1, FileType='Ascii')

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.4999999850988388, 0.4999999850988388, 1.5029607476163849]
renderView1.CameraFocalPoint = [0.4999999850988388, 0.4999999850988388, 0.4999999850988388]
renderView1.CameraParallelScale = 0.2595853468300873

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
