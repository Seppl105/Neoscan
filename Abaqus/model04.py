# coding=utf-8
from abaqus import *
from abaqusConstants import *
import regionToolset

# the positive x axis coincides with the r-direction in the modeled plane; the revolution wil be arounnd the y-axis

roundToDigit = 9 # if a value is rounded it will be to the "roundToDigit" digits

windingDivisor = 60       # Es wird für 600/windingDivisor Windungen gerechnet
totalWindings = int(600/windingDivisor)

#   Materialparameter des Leiterbands
t_con = 0.12 * 10**(-3) # [m] Dicke des Bandleiters
t_cow = 0.23 * 10**(-3) # [m] Dicke des Cowindings
t_ins = 0.01 * 10**(-3) # [m] Dicke der Isolation
t_total = t_con + t_cow + t_ins             # [m] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
materialWidths = [t_con, t_cow, t_ins]

density_con = 1.0
density_cow = 1.0
density_ins = 1.0

E_con = 500 * 10**9#280 * 10**9 # [Pa] E-Modul des Conductors
E_cow = 450 * 10**9#300 * 10**9 # [Pa] E-Modul des Cowindingsab
E_ins = 400 * 10**9#200 * 10**9 # [Pa] E-Modul der Insulation
materialEs = [E_con, E_cow, E_ins]

ny_con = 0.35 # Possion's Ratio Conductor
ny_cow = 0.3  # Possion's Ratio Cowinding
ny_ins = 0.4  # Possion's Ratio Insulation
materialNys = [ny_con, ny_cow, ny_ins]


r_i = 0.430               # [m] innerer Radius
r_a = 0.646               # [m] äußerer Radius ohne Confinement
if r_a != (r_i + t_total * ( ((r_a - r_i) / t_total) / windingDivisor)):
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n r_a und t Werte mit der Anzahl der Windungen ergeben nicht dev vorgegegebenen äußeren Radius\n vorgegebener Wert für r_a: ", r_a, "\n mit Dicken berechneter Wert round((r_i + t_total * ( ((r_a - r_i) / t_total) / windingDivisor), roundToDigit): ", round(r_i + t_total * ( ((r_a - r_i) / t_total) / windingDivisor), roundToDigit), "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
r_a = round(r_i + t_total * ( ((r_a - r_i) / t_total) / windingDivisor), roundToDigit)
yCoorStart = 0.0           # [m] x Coordinate am Anfang des ersten Pancakes
yCoorEnd = 0.01            # [m] x Coordinate am Ende des ersten Pancakes
yCoorMid = (yCoorEnd + yCoorStart)/2

sheetSizeCalculated = r_i + totalWindings*t_total + 5*t_total


# Create a new model
model = mdb.Model(name='Model-Coil')

# Define materials in Abaqus
materialCon = model.Material(name='material-Con')
materialCon.Density(table=((density_con, ), ))
materialCon.Elastic(table=((E_con, ny_con), ))

materialCow = model.Material(name='material-Cow')
materialCow.Density(table=((density_cow, ), ))
materialCow.Elastic(table=((E_cow, ny_cow), ))

materialIns = model.Material(name='material-Ins')
materialIns.Density(table=((density_ins, ), ))
materialIns.Elastic(table=((E_ins, ny_ins), ))

# Define sections for the material
model.HomogeneousSolidSection(name='Section-Con', material='Material-Con', thickness=None)
model.HomogeneousSolidSection(name='Section-Cow', material='Material-Cow', thickness=None)
model.HomogeneousSolidSection(name='Section-Ins', material='Material-Ins', thickness=None)
sectionNames = ['Section-Con', 'Section-Cow', 'Section-Ins'] # this array is used to assign the sectionsm, therfore the order is important and has to match materialWidths


# Create a new sketch
sketchOutline = model.ConstrainedSketch(name='Sketch-Coil', sheetSize=sheetSizeCalculated)
g, v, d, c = sketchOutline.geometry, sketchOutline.vertices, sketchOutline.dimensions, sketchOutline.constraints

sketchOutline.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0)) 
sketchOutline.FixedConstraint(entity=g[2])

#sketchOutline = sketchOutline.rectangle(point1=(r_i, yCoorStart), point2=(r_a, yCoorEnd))
sketchOutline.rectangle(point1=(r_i, yCoorStart), point2=(r_a, yCoorEnd))
sketchOutline.unsetPrimaryObject()
# Create the Part based on the sketch
part = model.Part(name='Part-Coil', dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY) ####################### TWO_D_PLANAR zu axial wechseln glaube ich
part.BaseShell(sketch=sketchOutline)

sectionPoints = [[] for element in sectionNames] # a list containing one list for each section which each list containing one point within each corresponding face

outerR = r_i    # outer radius of the partioning (rises which each added layer)
for winding in range(totalWindings):
    for i,width in enumerate(materialWidths):
        if outerR + width != r_a:
            # patition the face above the most recently created layer
            faceToPartition = part.faces.findAt(((outerR + width/2, yCoorMid, 0.0),)) # 3D point needed?
            
            #print("Partitioning at face: {faceToPartition}, with points: {(outerR + width, yCoorStart), (outerR + width, yCoorEnd)}")
            part.PartitionFaceByShortestPath(faces=faceToPartition, point1= (outerR + width, yCoorStart, 0.0), point2= (outerR + width, yCoorEnd, 0.0)) #2D points needed?
            
            #part.PartitionFaceByShortestPath(faces=faceToPartition, point1=(outerR + width, yCoorStart), point2=(outerR + width, yCoorEnd))  # 2D points needed?

            # Debugging print statements
            print("Outer Radius: ", outerR, ", Width: ", width, ", Material Index: ", i)
            print("Partitioning points: (", outerR + width, ", ", yCoorStart,"), (", outerR + width, ", ", yCoorEnd, ")")
            
            #faceOfNewLayer = part.faces.findAt(((outerR + width/2, yCoorMid, 0.0)))
            #faceOfNewLayerSequence = part.faces.sequenceFromLabels(labels=[faceOfNewLayer.index])
            #regionOfNewLayer = regionToolset.Region(faces=faceOfNewLayerSequence)
            #part.SectionAssignment(region=regionOfNewLayer, sectionName=sectionNames[i])
        # assign a point within the lastly partioned layer to the corresponding list (used to assign sections to the faces) 
        sectionPoints[i].append((outerR + width/2, yCoorMid, 0.0))
        outerR = round(outerR + width, roundToDigit)


# assign the faces to the corresponding sections according to sectionPoints
for i,sectionName in enumerate(sectionNames):
    tuplePoints = tuple([(point,) for point in sectionPoints[i]]) # creating a tuple of tuple ==> tuplPoints => (((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),))
    facesOfSectionName = part.faces.findAt(*tuplePoints)  # ==> the argument must be tuple of the shape ((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),) ==> hence we must remove the outer tuple ==> hence, we passed the argument with '*' sign which takes care of this
    
    #faceOfNewLayerSequence = part.faces.sequenceFromLabels(labels=[faceOfNewLayer.index])
    regionOfSectionName= regionToolset.Region(faces=facesOfSectionName)
    part.SectionAssignment(region=regionOfSectionName, sectionName=sectionName)


# Save the model
mdb.saveAs(pathName='Model-Coil')



# from abaqus import *
# from abaqusConstants import *
# import regionToolset

# # Create a model
# model = mdb.Model(name='Model-1')

# # Define materials
# material1 = model.Material(name='Material-1')
# material1.Density(table=((1.0, ), ))
# material1.Elastic(table=((210000.0, 0.3), ))

# material2 = model.Material(name='Material-2')
# material2.Density(table=((2.0, ), ))
# material2.Elastic(table=((100000.0, 0.25), ))

# material3 = model.Material(name='Material-3')
# material3.Density(table=((3.0, ), ))
# material3.Elastic(table=((150000.0, 0.35), ))

# # Create a sketch for the base feature
# s = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
# s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
# s.FixedConstraint(entity=s.geometry[2])
# s.rectangle(point1=(0.0, 0.0), point2=(10.0, 10.0))

# # Create the part
# part = model.Part(name='Part-1', dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY)
# part.BaseShell(sketch=s)

# # Number of layers
# num_layers = 6
# layer_thickness = 10.0 / num_layers  # Assuming the height of the part is 10.0

# # Partition the part into layers
# for i in range(1, num_layers):
#     y_coord = i * layer_thickness
    
#     # Find the face at the y-coordinate
#     face = part.faces.findAt(((5.0, y_coord - layer_thickness / 2, 0.0),))
    
#     # Partition the face
#     part.PartitionFaceByShortestPath(faces=face, point1=(0.0, y_coord), point2=(10.0, y_coord))

# # Create sections for the materials
# model.HomogeneousSolidSection(name='Section-1', material='Material-1', thickness=None)
# model.HomogeneousSolidSection(name='Section-2', material='Material-2', thickness=None)
# model.HomogeneousSolidSection(name='Section-3', material='Material-3', thickness=None)

# # Assign materials to every third layer
# for i in range(num_layers):
#     faces = part.faces.findAt(((5.0, (i + 0.5) * layer_thickness, 0.0), ))
#     region = regionToolset.Region(faces=faces)
#     if i % 3 == 0:
#         part.SectionAssignment(region=region, sectionName='Section-1')
#     elif i % 3 == 1:
#         part.SectionAssignment(region=region, sectionName='Section-2')
#     else:
#         part.SectionAssignment(region=region, sectionName='Section-3')

# # Save the model
# mdb.saveAs(pathName='axisymmetric_model_with_layers')


# # Draw a rectangle
# sketch.rectangle(point1=(0.0, 0.0), point2=(10.0, 5.0))

# # Create a new part based on the sketch
# part = model.Part(name='Part-1', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
# part.BaseShell(sketch=sketch)

# # Optional: Save the model to a file
# mdb.saveAs(pathName='model_with_sketch')