# coding=utf-8
from abaqus import *
from abaqusConstants import *
import regionToolset

#####################C:\Users\SGnir\OneDrive\TU\Bachelorsemester06\Kernspintomographie (MRT)\GitHub-Neoscan\Neoscan\Abaqus\model06.py
    

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

#   Materialparameter des Confinements
t_confinement = 0.0005       # [m] Dicke des Confinements (für "-1" wird ohne Confinement gerechnet)
E_confinement = E_con        # [Pa] E-Modul des Confinements
ny_confinement = ny_con      # Possions's Ratio des Confinements
density_confinement = 1.0


r_i = 0.430               # [m] innerer Radius
r_a = 0.646               # [m] äußerer Radius ohne Confinement
if r_a != (r_i + t_total * ( ((r_a - r_i) / t_total) / windingDivisor)):
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n r_a und t Werte mit der Anzahl der Windungen ergeben nicht dev vorgegegebenen äußeren Radius\n vorgegebener Wert für r_a: ", r_a, "\n mit Dicken berechneter Wert round((r_i + t_total * ( ((r_a - r_i) / t_total) / windingDivisor), roundToDigit): ", round(r_i + t_total * ( ((r_a - r_i) / t_total) / windingDivisor), roundToDigit), "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
r_a = round(r_i + t_total * ( ((r_a - r_i) / t_total) / windingDivisor), roundToDigit)
yCoorStart = 0.0           # [m] y Coordinate am Anfang des ersten Pancakes
yCoorEnd = 0.004            # [m] y Coordinate am Ende des ersten Pancakes
yCoorMid = (yCoorEnd + yCoorStart)/2

j = 114.2 * 10**6                        # [A/m^2] Stromdichte
b_za = -2.5                              # [T] magnetische Flussdichte außen
b_zi = 14                                # [T] magnetische Flussdichte innen
b_0 = b_za - (b_za-b_zi)/(r_a-r_i) * r_a # [T] absolutes Glied der Geradengleichung für B(r)

# "sheetSizeCalculated" is used to set the estimated size of the part before drawing
if t_confinement == -1:
    sheetSizeCalculated = r_i + totalWindings*t_total + 5*t_total
else:
    sheetSizeCalculated = r_i + totalWindings*t_total + 5*t_total + t_confinement


# Create a new model
model = mdb.Model(name='Model-Coil')

# Define materials in Abaqus
materialCon = model.Material(name='Material-Con')
materialCon.Density(table=((density_con, ), ))
materialCon.Elastic(table=((E_con, ny_con), ))

materialCow = model.Material(name='Material-Cow')
materialCow.Density(table=((density_cow, ), ))
materialCow.Elastic(table=((E_cow, ny_cow), ))

materialIns = model.Material(name='Material-Ins')
materialIns.Density(table=((density_ins, ), ))
materialIns.Elastic(table=((E_ins, ny_ins), ))


# Define sections for the material
model.HomogeneousSolidSection(name='Section-Con', material='Material-Con', thickness=None)
model.HomogeneousSolidSection(name='Section-Cow', material='Material-Cow', thickness=None)
model.HomogeneousSolidSection(name='Section-Ins', material='Material-Ins', thickness=None)
sectionNames = ['Section-Con', 'Section-Cow', 'Section-Ins'] # this array is used to assign the sectionsm, therfore the order is important and has to match to the widhts in "materialWidths"


# Create the sketch for one Pancake
sketchPancake = model.ConstrainedSketch(name='Sketch-Pancake', sheetSize=sheetSizeCalculated)
g, v, d, c = sketchPancake.geometry, sketchPancake.vertices, sketchPancake.dimensions, sketchPancake.constraints

###################### defining the axis of revolution?????????????????????????????????????????
#sketchPancake.sketchOptions.setValues(viewStyle=AXISYM) #######???
sketchPancake.ConstructionLine(point1=(0.0, - sheetSizeCalculated/2), point2=(0.0, sheetSizeCalculated/2)) 
sketchPancake.FixedConstraint(entity=g[2]) ################# 


sketchPancake.rectangle(point1=(r_i, yCoorStart), point2=(r_a, yCoorEnd))
#sketchPancake.unsetPrimaryObject()

# Create Part-Pancake based on the sketch
partPancake = model.Part(name='Part-Pancake', dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY)
partPancake.BaseShell(sketch=sketchPancake)

sectionPoints = [[] for element in sectionNames] # a list containing one list for each section with each list containing one point within each corresponding face
                                                 # ("sectionPoints" is later used to assign the sections)
                                                 
outerR = r_i    # outer radius of the partioning (rises which each added layer)
for winding in range(totalWindings):
    for i,width in enumerate(materialWidths):
        if outerR + width != r_a:
            # patition the face above the most recently created layer
            faceToPartition = partPancake.faces.findAt(((outerR + width/2, yCoorMid, 0.0),)) # 3D point needed?
            
            #print("Partitioning at face: {faceToPartition}, with points: {(outerR + width, yCoorStart), (outerR + width, yCoorEnd)}")
            partPancake.PartitionFaceByShortestPath(faces=faceToPartition, point1= (outerR + width, yCoorStart, 0.0), point2= (outerR + width, yCoorEnd, 0.0)) #2D points needed?
            
            #part.PartitionFaceByShortestPath(faces=faceToPartition, point1=(outerR + width, yCoorStart), point2=(outerR + width, yCoorEnd))  # 2D points needed?

            # Debugging print statements
            print("Outer Radius: ", outerR, ", Width: ", width, ", Material Index: ", i)
            print("Partitioning points: (", outerR + width, ", ", yCoorStart,"), (", outerR + width, ", ", yCoorEnd, ")")
            
            #faceOfNewLayer = partPancake.faces.findAt(((outerR + width/2, yCoorMid, 0.0)))
            #faceOfNewLayerSequence = partPancake.faces.sequenceFromLabels(labels=[faceOfNewLayer.index])
            #regionOfNewLayer = regionToolset.Region(faces=faceOfNewLayerSequence)
            #partPancake.SectionAssignment(region=regionOfNewLayer, sectionName=sectionNames[i])
        # assign a point within the lastly partioned layer to the corresponding list (used to assign sections to the faces) 
        sectionPoints[i].append((outerR + width/2, yCoorMid, 0.0))
        outerR = round(outerR + width, roundToDigit)

print("outerR equal to assigned or calculated r_a: ", outerR == r_a, " outerR: ", outerR, " r_a: ", r_a)

# assign the faces to the corresponding sections according to sectionPoints
for i,sectionName in enumerate(sectionNames):
    tuplePoints = tuple([(point,) for point in sectionPoints[i]]) # creating a tuple of tuple ==> tuplPoints => (((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),))
    facesOfSectionName = partPancake.faces.findAt(*tuplePoints)  # ==> the argument must be tuple of the shape ((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),) ==> hence we must remove the outer tuple ==> hence, we passed the argument with '*' sign which takes care of this
    
    #faceOfNewLayerSequence = partPancake.faces.sequenceFromLabels(labels=[faceOfNewLayer.index])
    regionOfSectionName= regionToolset.Region(faces=facesOfSectionName)
    partPancake.SectionAssignment(region=regionOfSectionName, sectionName=sectionName)
    

# add Part-Pancake to assembly
assembly = mdb.models['Model-Coil'].rootAssembly
assembly.DatumCsysByThreePoints(coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0), 
    point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, -1.0))

partPancake = mdb.models['Model-Coil'].parts['Part-Pancake'] # may be redundant

assembly.Instance(name='Part-Pancake-AssemblyInstance', part=partPancake, dependent=ON)

if t_confinement != -1:
    # create Material for the confinement
    materialConfinement = model.Material(name='Material-Confinement')
    materialConfinement.Density(table=((density_confinement, ), ))
    materialConfinement.Elastic(table=((E_confinement, ny_confinement), ))
    # creating Section-Confinement with Material-Confinement
    model.HomogeneousSolidSection(name='Section-Confinement', material='Material-Confinement', thickness=None)
    
    # create Sketch-Confinement
    sketchConfinement = mdb.models['Model-Coil'].ConstrainedSketch(name='Sketch-Confinement', sheetSize=sheetSizeCalculated)
    g, v, d, c = sketchConfinement.geometry, sketchConfinement.vertices, sketchConfinement.dimensions, sketchConfinement.constraints
    sketchConfinement.sketchOptions.setValues(viewStyle=AXISYM)
    sketchConfinement.setPrimaryObject(option=STANDALONE) # may be needless
    sketchConfinement.ConstructionLine(point1=(0.0, - sheetSizeCalculated/2), point2=(0.0, sheetSizeCalculated/2))
    sketchConfinement.FixedConstraint(entity=g[2])
    yOffSetConfinement = 0#20*t_confinement   # this is used to create a visible offset
    xOffSetConfinement = 20*t_confinement   # this is used to create a visible offset
    sketchConfinement.rectangle(point1=(r_a + xOffSetConfinement, yCoorStart + yOffSetConfinement), point2=(r_a + xOffSetConfinement + t_confinement, yCoorEnd + yOffSetConfinement))
    
    # create Part-Confinment
    partConfinement = mdb.models['Model-Coil'].Part(name='Part-Confinement', dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY)
    partConfinement.BaseShell(sketch=sketchConfinement)
    #sketchConfinement.unsetPrimaryObject()  # may be needless
    #p = mdb.models['Model-Coil'].parts['Part-Confinement']
    #session.viewports['Viewport: 1'].setValues(displayedObject=p)
    #del mdb.models['Model-Coil'].sketches['__profile__']
    ###assembly = mdb.models['Model-Coil'].rootAssembly
    #session.viewports['Viewport: 1'].setValues(displayedObject=a)
    #session.viewports['Viewport: 1'].view.setValues(nearPlane=0.857601, farPlane=0.877863, width=0.0584114, height=0.0239605, viewOffsetX=0.188462, viewOffsetY=0.000243931)

    # adding Part-Confinement to Section-Confinement which is linked to Material-Confinement
    faceConfinement = partConfinement.faces.findAt(((r_a + t_confinement/2 + xOffSetConfinement, yCoorMid + yOffSetConfinement, 0.0),))
    region = partConfinement.Set(faces=faceConfinement, name='Set-Confinement')
    partConfinement = mdb.models['Model-Coil'].parts['Part-Confinement'] # may be redundant
    partConfinement.SectionAssignment(region=region, sectionName='Section-Confinement', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    
    # add Part-Confinement to the assembly
    assembly = mdb.models['Model-Coil'].rootAssembly # may be redundant
    partConfinement = mdb.models['Model-Coil'].parts['Part-Confinement']
    assembly.Instance(name='Part-Confinement-AssemblyInstance', part=partConfinement, dependent=ON)
    
    # set a constraint which coincides the outer edge of the pancake to the confinement
    #edgeMovable = mdb.models['Model-Coil'].parts['Part-Confinement'].edges.getByBoundingBox()
    #edge2 = assembly.instances['Part-Pancake-AssemblyInstance'].edges
    
    
    print("edge of Part-Confinement: ", mdb.models['Model-Coil'].parts['Part-Confinement'].edges.findAt(((r_a + xOffSetConfinement, yCoorMid + yOffSetConfinement, 0))))
    #print("edge: ", mdb.models['Model-Coil'].parts['Part-Confinement'].edges.findAt(((r_a + xOffSetConfinement, yCoorMid + yOffSetConfinement, 0)))[0])
    print("edge of Part-Pancake: ", mdb.models['Model-Coil'].parts['Part-Pancake'].edges.findAt(((r_a, yCoorMid, 0))))
    
    #edge1 = assembly.instances['Part-Confinement-AssemblyInstance'].edges
    #print("assembly.instances['Part-Confinement-AssemblyInstance'].InterestingPoint(edge=edge1[3], rule=MIDDLE), : ", assembly.instances['Part-Confinement-AssemblyInstance'].InterestingPoint(edge=edge1[3], rule=MIDDLE))
    
    #movabelEdge = mdb.models['Model-Coil'].parts['Part-Confinement'].edges.findAt(((r_a + xOffSetConfinement, yCoorMid + yOffSetConfinement, 0)))
    
    movabelEdge = assembly.instances['Part-Confinement-AssemblyInstance'].edges.findAt(((r_a + xOffSetConfinement, yCoorMid + yOffSetConfinement, 0)))
    fixedEdge = assembly.instances['Part-Pancake-AssemblyInstance'].edges.findAt(((r_a, yCoorMid, 0)))
    
    print("assembly.instances['Part-Confinement-AssemblyInstance']: ", assembly.instances['Part-Confinement-AssemblyInstance'].InterestingPoint(edge=movabelEdge, rule=MIDDLE))
    #print("assembly.instances['Part-Pancake-AssemblyInstance']: ", assembly.instances['Part-Pancake-AssemblyInstance'].InterestingPoint(edge=mdb.models['Model-Coil'].parts['Part-Pancake'].edges.findAt(((r_a, yCoorMid, 0))), rule=MIDDLE))
    #print("edge of Part-Pancake: ", mdb.models['Model-Coil'].parts['Part-Pancake-1'].edges.findAt(((r_a, yCoorMid, 0)))) # stops the flow with an error
    assembly.CoincidentPoint(
        movablePoint=assembly.instances['Part-Confinement-AssemblyInstance'].InterestingPoint(edge=movabelEdge, rule=MIDDLE), 
        fixedPoint=assembly.instances['Part-Pancake-AssemblyInstance'].InterestingPoint(edge=fixedEdge, rule=MIDDLE))
    print("worked")
    #print("edge of Part-Pancake: ", mdb.models['Model-Coil'].parts['Part-Pancake-1'].edges.findAt(((r_a, yCoorMid, 0)))) # stops the flow with an error
    # edge1 = assembly.instances['Part-Confinement-AssemblyInstance'].edges
    # edge2 = assembly.instances['Part-Pancake-AssemblyInstance'].edges
    # assembly.CoincidentPoint(
    #     movablePoint=assembly.instances['Part-Confinement-AssemblyInstance'].InterestingPoint(
    #     edge=edge1[3], rule=MIDDLE), 
    #     fixedPoint=assembly.instances['Part-Pancake-AssemblyInstance'].InterestingPoint(edge=edge2[89], 
    #     rule=MIDDLE))
    
    # create a contact(/intersection) property
    mdb.models['Model-Coil'].ContactProperty('IntProp-1')
    # set the coefficient of friction to infinity (by using "ROUGH")
    mdb.models['Model-Coil'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=ROUGH)
    # prevent overlaping by using "HARD" and allow Seperation
    mdb.models['Model-Coil'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    # assign contact propertys to all intersections (reminder: intersections are edges between parts)
    mdb.models['Model-Coil'].ContactStd(name='Int-1', createStepName='Initial')
    mdb.models['Model-Coil'].interactions['Int-1'].includedPairs.setValuesInStep(
        stepName='Initial', useAllstar=ON)
    mdb.models['Model-Coil'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
        stepName='Initial', assignments=((GLOBAL, SELF, 'IntProp-1'), ))




# Create StaticStep named Step-1-BodyForce
mdb.models['Model-Coil'].StaticStep(name='Step-1-BodyForce', previous='Initial', description='Description-Step-1-Body-Force')


# Create Load (body force)
expressionStringBodyForce = str(j) + ' * ( ((' + str(b_za) + ') - (' + str(b_zi) + '))/(' + str(r_a) + '-' + str(r_i) + ') * X + (' + str(b_za) + ') - ((' + str(b_za) + ') - (' + str(b_zi) + '))/(' + str(r_a) + '-' + str(r_i) + ') * ' + str(r_a) + ')' # '114.2*pow ( 10,6) *( ((-2.5)- (14))/(0.646-0.430) * X + (-2.5) - ((-2.5) - (14))/(0.646-0.430)  * 0.646)'
print("Expression to calculate the body Force dependant on X: ", expressionStringBodyForce, "\n maybe helpful to compare to: ", '114.2*pow ( 10,6) *( ((-2.5)- (14))/(0.646-0.430) * X + (-2.5) - ((-2.5) - (14))/(0.646-0.430)  * 0.646)')
mdb.models['Model-Coil'].ExpressionField(name='AnalyticalField-1-BodyForce', localCsys=None, description='Description-AnalyticalField-1-BodyForce', expression=expressionStringBodyForce)
#mdb.models['Model-Coil'].ExpressionField(name='AnalyticalField-1-BodyForce', localCsys=None, description='Description-AnalyticalField-1-BodyForce', 
#                            expression='114.2*pow ( 10,6) * ((-2.5)- (14))/(0.646-0.430) * X + (-2.5) - ((-2.5) - (14))/(0.646-0.430)  * 0.646)')
#mdb.models['Model-Coil'].TabularAmplitude(name='Amp-1', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 1.0), (1.0, 1.0)))############################### step time for static step?
assembly = mdb.models['Model-Coil'].rootAssembly # may be redundant
facesPancake = assembly.instances['Part-Pancake-AssemblyInstance'].faces

facesPancake = facesPancake.getSequenceFromMask(mask=('[#3fffffff ]', ), )
regionFacesOfPancake = assembly.Set(faces=facesPancake, name='Set-FacesPancake')
mdb.models['Model-Coil'].BodyForce(name='Load-1-BodyForce', createStepName='Step-1-BodyForce', region=regionFacesOfPancake, comp1=1.0, distributionType=FIELD, field='AnalyticalField-1-BodyForce')



# pointsInPancake = []
# print(sectionPoints)
# for listOfPoints in sectionPoints:
#     print(listOfPoints)
#     pointsInPancake.extend(listOfPoints)
# print("pointsInPancake: ", pointsInPancake )
# #pointsInPancake = [pointsInPancake.extend(listOfPoints) for listOfPoints in sectionPoints]
# print(pointsInPancake)
# tuplePointsInPancake = tuple([(point,) for point in pointsInPancake]) # creating a tuple of tuple ==> tuplPoints => (((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),))
# print(tuplePointsInPancake)
# facesOfPancake = partPancake.faces.findAt(*tuplePointsInPancake)  # ==> the argument must be tuple of the shape ((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),) ==> hence we must remove the outer tuple ==> hence, we passed the argument with '*' sign which takes care of this
# #facesOfPancake = f1.getSequenceFromMask(mask=('[#3fffffff ]', ), )
# #regionFacesOfPancake = assembly.Set(faces=facesOfPancake, name='Set-1-BodyForce')
# print("ttttttttttttttttttttttt")
# regionFacesOfPancake = regionToolset.Region(faces=facesOfPancake)
# print("test")
# #mdb.models['Model-Coil'].BodyForce(name='Load-1-BodyForce', createStepName='Step-1-BodyForce', region=regionFacesOfPancake, comp1=1.0, distributionType=FIELD, field='AnalyticalField-1-BodyForce')



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