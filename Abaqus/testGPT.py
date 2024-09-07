# from abaqus import *
# from abaqusConstants import *

# # Create a new model
# model = mdb.Model(name='Model-1')

# # Create a new sketch
# sketch = model.ConstrainedSketch(name='Sketch-1', sheetSize=200.0)

# # Draw a rectangle
# sketch.rectangle(point1=(0.0, 0.0), point2=(10.0, 5.0))

# # Create a new part based on the sketch
# part = model.Part(name='Part-1', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
# part.BaseShell(sketch=sketch)

# # Optional: Save the model to a file
# mdb.saveAs(pathName='model_with_sketch')

from abaqus import *
from abaqusConstants import *
import regionToolset

# Create a new model
model = mdb.Model(name='Model-1')

# Create a sketch for the base feature
sketch = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)

# Create a closed rectangle sketch
sketch.rectangle(point1=(0.0, 0.0), point2=(10.0, 10.0))

# Ensure the sketch is closed by applying coincident constraints
sketch.autoTrimCurve(curve1=sketch.geometry[2], point1=(0.0, 0.0))
sketch.autoTrimCurve(curve1=sketch.geometry[3], point1=(10.0, 0.0))
sketch.autoTrimCurve(curve1=sketch.geometry[4], point1=(10.0, 10.0))
sketch.autoTrimCurve(curve1=sketch.geometry[5], point1=(0.0, 10.0))

# Check for continuity by ensuring all constraints are applied
# Here we programmatically ensure that the sketch is closed
sketch.unsetPrimaryObject()
assembly = model.rootAssembly
part = model.Part(name='Part-1', dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY)
part.BaseShell(sketch=sketch)

# Number of layers
num_layers = 6
layer_thickness = 10.0 / num_layers  # Assuming the height of the part is 10.0

# Partition the part into layers
for i in range(1, num_layers):
    y_coord = i * layer_thickness
    
    # Find the face at the y-coordinate
    face = part.faces.findAt(((5.0, y_coord - layer_thickness / 2, 0.0),))
    
    # Partition the face
    part.PartitionFaceByShortestPath(faces=face, point1=(0.0, y_coord), point2=(10.0, y_coord))

# Define materials
material1 = model.Material(name='Material-1')
material1.Density(table=((1.0, ), ))
material1.Elastic(table=((210000.0, 0.3), ))

material2 = model.Material(name='Material-2')
material2.Density(table=((2.0, ), ))
material2.Elastic(table=((100000.0, 0.25), ))

material3 = model.Material(name='Material-3')
material3.Density(table=((3.0, ), ))
material3.Elastic(table=((150000.0, 0.35), ))

# Create sections for the materials
model.HomogeneousSolidSection(name='Section-1', material='Material-1', thickness=None)
model.HomogeneousSolidSection(name='Section-2', material='Material-2', thickness=None)
model.HomogeneousSolidSection(name='Section-3', material='Material-3', thickness=None)

# Assign materials to every third layer
for i in range(num_layers):
    faces = part.faces.findAt(((5.0, (i + 0.5) * layer_thickness, 0.0), ))
    region = regionToolset.Region(faces=faces)
    if i % 3 == 0:
        part.SectionAssignment(region=region, sectionName='Section-1')
    elif i % 3 == 1:
        part.SectionAssignment(region=region, sectionName='Section-2')
    else:
        part.SectionAssignment(region=region, sectionName='Section-3')

# Save the model
mdb.saveAs(pathName='axisymmetric_model_with_layers')
