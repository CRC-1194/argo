#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
from numpy import random

# create a new 'Cylinder'
cylinder1 = Cylinder()

# Properties modified on cylinder1
cylinder1.Resolution = 512
cylinder1.Height = 0.1

# Random cylinder positioning 
cylinder1.Radius = 0.15

# Manually perturbed position,
cylinder1.Center = [0, 0, 0]

cylinder1.Capping = 0

# create a new 'Transform'
transform1 = Transform(Input=cylinder1)
transform1.Transform = 'Transform'

# Properties modified on transform1.Transform
transform1.Transform.Rotate = [90.0, 0.0, 0.0]

# create a new 'Transform'
transform2 = Transform(Input=transform1)
transform2.Transform = 'Transform'

translation = [0.5,  0.5,  0.5 * cylinder1.Height]

# Properties modified on transform2.Transform
transform2.Transform.Translate = translation

# create a new 'Triangulate'
triangulate1 = Triangulate(Input=transform2)

# Save cylinder data
SaveData('./cylinder.stl', proxy=triangulate1, FileType='Ascii')

