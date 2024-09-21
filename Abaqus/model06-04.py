# coding=utf-8
from abaqus import *
from abaqusConstants import *
import regionToolset
import displayGroupOdbToolset as dgo
import displayGroupMdbToolset as dgm
#####################C:\Users\SGnir\OneDrive\TU\Bachelorsemester06\Kernspintomographie (MRT)\GitHub-Neoscan\Neoscan\Abaqus\model06.py
    

# the positive x axis coincides with the r-direction in the modeled plane; the revolution wil be arounnd the y-axis
loadOnAll = False
roundToDigit = 9 # if a value is rounded it will be to the "roundToDigit" digits

# windingDivisor = 60       # Es wird für 600/windingDivisor Windungen gerechnet
# totalWindings = int(600/windingDivisor)
totalWindings = 5 #600
windingDivisor = 600/totalWindings

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
yCoorEnd = 0.004#0.00036 #0.004            # [m] y Coordinate am Ende des ersten Pancakes
yCoorMid = (yCoorEnd + yCoorStart)/2
yReal = 0.004

j = 250/(yReal*t_con)    #114.2 * 10**6                        # [A/m^2] Stromdichte
b_za = -2.5                              # [T] magnetische Flussdichte außen
b_zi = 14                                # [T] magnetische Flussdichte innen
b_0 = b_za - (b_za-b_zi)/(r_a-r_i) * r_a # [T] absolutes Glied der Geradengleichung für B(r)

# "sheetSizeCalculated" is used to set the estimated size of the part before drawing
# if t_confinement == -1:
#     sheetSizeCalculated = r_i + totalWindings*t_total + 5*t_total
# else:
#     sheetSizeCalculated = r_i + totalWindings*t_total + 5*t_total + t_confinement
#sheetSizeCalculated = t_total*10 + t_confinement

# Create a new model
model = mdb.Model(name='Model-1')

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

parts = ["Con", "Cow", "Ins"]

def createRing(name, pointOne, pointTwo):
    # points have to be tuples e.g. (0.0, 2.0)
    sketchRing = model.ConstrainedSketch(name="Sketch-" + name, sheetSize=max(abs(pointOne[0] - pointTwo[0]), abs(pointOne[1] - pointTwo[1])))
    g, v, d, c = sketchRing.geometry, sketchRing.vertices, sketchRing.dimensions, sketchRing.constraints

    
    ###################### defining the axis of revolution?????????????????????????????????????????
    sketchRing.ConstructionLine(point1=(0.0, pointOne[1]), point2=(0.0, pointTwo[1])) 
    sketchRing.FixedConstraint(entity=g[2]) ################# 
    
    sketchRing.rectangle(point1=pointOne, point2=pointTwo)
    
    ring = model.Part(name="Part-" + name, dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY)
    ring.BaseShell(sketch=sketchRing)
    
    del model.sketches["Sketch-" + name] #free up memory space
    
    return ring

# put all that into the loop as wellll !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Create the first ring
part = createRing(name=parts[0] + "1", pointOne=(r_i, yCoorStart), pointTwo=(r_i + materialWidths[0], yCoorEnd))
# assign section
faceOfPart = part.faces.findAt(((r_i + materialWidths[0]/2, yCoorMid, 0.0), ))  # the argument must be tuple of the shape ((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),)
regionOfSectionName= regionToolset.Region(faces=faceOfPart)
part.SectionAssignment(region=regionOfSectionName, sectionName=sectionNames[0])

# add the first ring to the assembly
assembly = mdb.models['Model-1'].rootAssembly
assembly.Instance(name="AssemblyInstance-" + parts[0] + "1", part=part, dependent=OFF)

# create "instanceNames" to Mesh all parts in the end
instanceNames = ["AssemblyInstance-" + parts[0] + "1"]

# sectionPoints = [[] for element in sectionNames] # a list containing one list for each section with each list containing one point within each corresponding face
#                                                  # ("sectionPoints" is later used to assign the sections)                                             

# Create all remaining parts, add them to the assembly and update sectionPoints
currentR = r_i    # outer radius of the partioning (rises which each added layer)
for winding in range(totalWindings):
    for i,width in enumerate(materialWidths):
        if winding == 0 and i == 0:
            # Used to assign the BC and the Node Sets on the Bottom later on
            edges = assembly.instances["AssemblyInstance-" + parts[0] + "1"].edges
            edgesBottom = edges.findAt(((r_i + materialWidths[0]/2, yCoorStart, 0.0), )) # all edges on the x-Axis will be stored within this variable
            edgesTop = edges.findAt(((r_i + materialWidths[0]/2, yCoorEnd, 0.0), )) ##################################################################################

            # Used to assign the Load later on
            faces = assembly.instances["AssemblyInstance-" + parts[0] + "1"].faces
            facesBodyForce = faces.findAt(((r_i + materialWidths[0]/2, yCoorMid, 0.0), ))

            # Used to assign the interaction properties later on
            edges = assembly.instances["AssemblyInstance-" + parts[0] + "1"].edges
            edgesMaster = edges.findAt(((r_i, yCoorMid, 0.0), )) # this face is just a dummy to create the variable, the face is not touching any other face!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            edgesSlave = edges.findAt(((r_i + materialWidths[0], yCoorMid, 0.0), ))

        else:
            # Create ring
            part = createRing(name=parts[i] + str(winding + 1), pointOne=(currentR, yCoorStart), pointTwo=(currentR + width, yCoorEnd))
            # Add ring to the assembly
            instanceName = "AssemblyInstance-" + parts[i] + str(winding + 1)
            assembly.Instance(name=instanceName, part=part, dependent=OFF)
            
            # add Instance Name to "instanceNames" to Mesh all Instances in the End
            instanceNames.append(instanceName)
            
            # assign section
            faceOfPart = part.faces.findAt(((currentR + width/2, yCoorMid, 0.0), ))  # the argument must be tuple of the shape ((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),)
            regionOfSectionName= regionToolset.Region(faces=faceOfPart)
            part.SectionAssignment(region=regionOfSectionName, sectionName=sectionNames[i])
            # old:
            # # assign a point within each instance to the corresponding list (used to assign sections) 
            # sectionPoints[i].append((currentR + width/2, yCoorMid, 0.0))
            
            # add edges to "edgesBottom" to assign the BC and the Node Sets later on
            edges = assembly.instances[instanceName].edges
            edgesBottom = edgesBottom + edges.findAt(((currentR + width/2, yCoorStart, 0.0),))
            edgesTop = edgesTop + edges.findAt(((currentR + width/2, yCoorEnd, 0.0),))####################################################################################
            #edgesBottom.append(edges.findAt(((currentR + width/2, yCoorStart, 0.0),)))
            
            
            # add faces to "facesBodyForce" to assign the Load later on
            if loadOnAll or i==0:
                faces = assembly.instances[instanceName].faces
                facesBodyForce = facesBodyForce + faces.findAt(((currentR + width/2, yCoorMid, 0.0),))
            
            # add edges to "edgesMaster" and "edgesSlave" to assign interaction conditions later on
            edges = assembly.instances[instanceName].edges
            if i == 0:
                edgesMaster = edgesMaster + edges.findAt(((currentR, yCoorMid, 0.0),))
                edgesSlave = edgesSlave + edges.findAt(((currentR + width, yCoorMid, 0.0),))
            elif i == 1:
                edgesMaster = edgesMaster + edges.findAt(((currentR, yCoorMid, 0.0),)) + edges.findAt(((currentR + width, yCoorMid, 0.0),))
            elif i == 2:
                edgesSlave = edgesSlave + edges.findAt(((currentR, yCoorMid, 0.0),)) + edges.findAt(((currentR + width, yCoorMid, 0.0),))

            # edges of the different materials for the mesh
            if i == 0:
                edgesMaster = edgesMaster + edges.findAt(((currentR, yCoorMid, 0.0),))
                edgesSlave = edgesSlave + edges.findAt(((currentR + width, yCoorMid, 0.0),))
            elif i == 1:
                edgesMaster = edgesMaster + edges.findAt(((currentR, yCoorMid, 0.0),)) + edges.findAt(((currentR + width, yCoorMid, 0.0),))
            elif i == 2:
                edgesSlave = edgesSlave + edges.findAt(((currentR, yCoorMid, 0.0),)) + edges.findAt(((currentR + width, yCoorMid, 0.0),))
                
        currentR = round(currentR + width, roundToDigit)
    # Used to assign the Nodesets later on
    if winding == 0:
        edges = assembly.instances["AssemblyInstance-" + parts[0] + str(winding+1)].edges
        edgesCon = edges.findAt(((currentR - t_ins - t_cow - t_con, yCoorMid, 0.0), ))
        edges = assembly.instances["AssemblyInstance-" + parts[1] + str(winding+1)].edges
        edgesCow = edges.findAt(((currentR - t_ins - t_cow, yCoorMid, 0.0), ))
        edges = assembly.instances["AssemblyInstance-" + parts[2] + str(winding+1)].edges
        edgesIns = edges.findAt(((currentR - t_ins, yCoorMid, 0.0), ))  
    else:
        edges = assembly.instances["AssemblyInstance-" + parts[0] + str(winding+1)].edges
        edgesCon = edgesCon + edges.findAt(((currentR - t_ins - t_cow - t_con, yCoorMid, 0.0), ))
        edges = assembly.instances["AssemblyInstance-" + parts[1] + str(winding+1)].edges
        edgesCow = edgesCow + edges.findAt(((currentR - t_ins - t_cow, yCoorMid, 0.0), ))
        edges = assembly.instances["AssemblyInstance-" + parts[2] + str(winding+1)].edges
        edgesIns = edges + edges.findAt(((currentR - t_ins, yCoorMid, 0.0), ))      
            

print("currentR equal to assigned or calculated r_a: ", currentR == r_a, " currentR: ", currentR, " r_a: ", r_a)
assembly.regenerate()

# Create StaticStep named Step-1-BodyForce
mdb.models['Model-1'].StaticStep(name='Step-1-BodyForce', previous='Initial', description='Description-Step-1-Body-Force')
mdb.models['Model-1'].steps['Step-1-BodyForce'].setValues(initialInc=0.05)



# create BC
regionBc = assembly.Set(edges=edgesBottom, name="Set-Bc")
mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Initial', 
                                     region=regionBc, u1=UNSET, u2=SET, ur3=UNSET, amplitude=UNSET, 
                                     distributionType=UNIFORM, fieldName='', localCsys=None)

# create Node Sets and Mesh
#assembly.seedEdgeByNumber(edges=edgesBottom, number=3, constraint=FINER)
meshStructure = 2
print("Mesh Structure: ", meshStructure)
meshWidthNumber = 4
meshHightToWidthRatio = 3
if meshStructure == 1:    
    assembly.seedEdgeBySize(edges=edgesBottom, size=3.33e-06, deviationFactor=0.1e-06, 
            minSizeFactor=0.1, constraint=FINER)
    assembly.seedEdgeBySize(edges=edgesSlave+edgesMaster, size=3.33e-06, deviationFactor=0.1e-06, 
            minSizeFactor=0.1, constraint=FINER)
elif meshStructure == 2:
    assembly.seedEdgeByNumber(edges=edgesBottom, number=meshWidthNumber, constraint=FINER)
    assembly.seedEdgeByNumber(edges=edgesTop, number=meshWidthNumber, constraint=FINER)
    assembly.seedEdgeBySize(edges=edgesCon, size=t_con/meshWidthNumber * meshHightToWidthRatio, deviationFactor=0.1e-06, minSizeFactor=0.1, constraint=FINER)
    assembly.seedEdgeBySize(edges=edgesCow, size=t_cow/meshWidthNumber * meshHightToWidthRatio,
            deviationFactor=0.1e-06, minSizeFactor=0.1, constraint=FINER)
    assembly.seedEdgeBySize(edges=edgesIns, size=t_ins/meshWidthNumber * meshHightToWidthRatio,
            deviationFactor=0.1e-06, minSizeFactor=0.1, constraint=FINER)
elif meshStructure == 3:
    assembly.seedEdgeByNumber(edges=edgesBottom, number=meshWidthNumber, constraint=FINER)
    assembly.seedEdgeByNumber(edges=edgesTop, number=meshWidthNumber, constraint=FINER)
    assembly.seedEdgeBySize(edges=edgesCon, size=t_con/meshWidthNumber,
                            deviationFactor=0.1e-06, minSizeFactor=0.1, constraint=FINER)
    assembly.seedEdgeBySize(edges=edgesCow, size=t_cow/meshWidthNumber,
            deviationFactor=0.1e-06, minSizeFactor=0.1, constraint=FINER)
    assembly.seedEdgeBySize(edges=edgesIns, size=t_ins/meshWidthNumber * meshHightToWidthRatio,
            deviationFactor=0.1e-06, minSizeFactor=0.1, constraint=FINER)
    
partInstances = tuple([assembly.instances[instanceName] for instanceName in instanceNames])
assembly.generateMesh(regions=partInstances)
    



# Create Load (body force)
expressionStringBodyForce = str(j) + ' * ( ((' + str(b_za) + ') - (' + str(b_zi) + '))/(' + str(r_a) + '-' + str(r_i) + ') * X + (' + str(b_za) + ') - ((' + str(b_za) + ') - (' + str(b_zi) + '))/(' + str(r_a) + '-' + str(r_i) + ') * ' + str(r_a) + ')' # '114.2*pow ( 10,6) *( ((-2.5)- (14))/(0.646-0.430) * X + (-2.5) - ((-2.5) - (14))/(0.646-0.430)  * 0.646)'
print("Expression to calculate the body Force dependant on X: ", expressionStringBodyForce, "\n maybe helpful to compare to: ", '114.2*pow ( 10,6) *( ((-2.5)- (14))/(0.646-0.430) * X + (-2.5) - ((-2.5) - (14))/(0.646-0.430)  * 0.646)')
mdb.models['Model-1'].ExpressionField(name='AnalyticalField-1-BodyForce', localCsys=None, description='Description-AnalyticalField-1-BodyForce', expression=expressionStringBodyForce)
regionBodyForce = assembly.Set(faces=facesBodyForce, name="Set-BodyForce")
#mdb.models['Model-Coil'].ExpressionField(name='AnalyticalField-1-BodyForce', localCsys=None, description='Description-AnalyticalField-1-BodyForce', 
#                            expression='114.2*pow ( 10,6) * ((-2.5)- (14))/(0.646-0.430) * X + (-2.5) - ((-2.5) - (14))/(0.646-0.430)  * 0.646)')
#mdb.models['Model-Coil'].TabularAmplitude(name='Amp-1', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 1.0), (1.0, 1.0)))############################### step time for static step?
#facesPancake = assembly.instances['Part-Pancake-AssemblyInstance'].faces
#facesPancake = facesPancake.getSequenceFromMask(mask=('[#3fffffff ]', ), )
#regionFacesOfPancake = assembly.Set(faces=facesPancake, name='Set-FacesPancake')
mdb.models['Model-1'].BodyForce(name='Load-1-BodyForce', createStepName='Step-1-BodyForce', region=regionBodyForce, comp1=1.0, distributionType=FIELD, field='AnalyticalField-1-BodyForce')

# create interactionProperties
mdb.models['Model-1'].ContactProperty('IntProp-1')
# mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
#     formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
#     pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
#     table=((0.9, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
#     fraction=0.005, elasticSlipStiffness=None)
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=FRICTIONLESS)
mdb.models['Model-1'].interactionProperties['IntProp-1'].tangentialBehavior.setValues(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
    table=((0.9, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
    fraction=0.005, elasticSlipStiffness=None)
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON, 
    constraintEnforcementMethod=DEFAULT)
# session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
# leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
# session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
# leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
# session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863083, 
#     farPlane=0.864223, width=0.00411057, height=0.00168617, 
#     viewOffsetX=0.214927, viewOffsetY=-0.00120241)
regionMaster = assembly.Surface(side1Edges=edgesMaster, name="m_Surf-Master")
regionSlave = assembly.Surface(side1Edges=edgesSlave, name="m_Surf-Slave")
# a = mdb.models['Model-1'].rootAssembly
# s1 = a.instances['AssemblyInstanceCow1'].edges
# side1Edges1 = s1.getSequenceFromMask(mask=('[#a ]', ), )
# s2 = a.instances['AssemblyInstanceCow2'].edges
# side1Edges2 = s2.getSequenceFromMask(mask=('[#a ]', ), )
# s3 = a.instances['AssemblyInstanceCow3'].edges
# side1Edges3 = s3.getSequenceFromMask(mask=('[#a ]', ), )
# s4 = a.instances['AssemblyInstanceCow4'].edges
# side1Edges4 = s4.getSequenceFromMask(mask=('[#a ]', ), )
# s5 = a.instances['AssemblyInstanceCow5'].edges
# side1Edges5 = s5.getSequenceFromMask(mask=('[#a ]', ), )
# region1=a.Surface(side1Edges=side1Edges1+side1Edges2+side1Edges3+side1Edges4+\
#     side1Edges5, name='m_Surf-1')
# a = mdb.models['Model-1'].rootAssembly
# s1 = a.instances['AssemblyInstance-Con1'].edges
# side1Edges1 = s1.getSequenceFromMask(mask=('[#2 ]', ), )
# s2 = a.instances['AssemblyInstanceIns1'].edges
# side1Edges2 = s2.getSequenceFromMask(mask=('[#8 ]', ), )
# s3 = a.instances['AssemblyInstanceCon2'].edges
# side1Edges3 = s3.getSequenceFromMask(mask=('[#2 ]', ), )
# s4 = a.instances['AssemblyInstanceIns2'].edges
# side1Edges4 = s4.getSequenceFromMask(mask=('[#8 ]', ), )
# s5 = a.instances['AssemblyInstanceCon3'].edges
# side1Edges5 = s5.getSequenceFromMask(mask=('[#2 ]', ), )
# s6 = a.instances['AssemblyInstanceIns3'].edges
# side1Edges6 = s6.getSequenceFromMask(mask=('[#8 ]', ), )
# s7 = a.instances['AssemblyInstanceCon4'].edges
# side1Edges7 = s7.getSequenceFromMask(mask=('[#2 ]', ), )
# s8 = a.instances['AssemblyInstanceIns4'].edges
# side1Edges8 = s8.getSequenceFromMask(mask=('[#8 ]', ), )
# s9 = a.instances['AssemblyInstanceCon5'].edges
# side1Edges9 = s9.getSequenceFromMask(mask=('[#2 ]', ), )
# s10 = a.instances['AssemblyInstanceIns5'].edges
# side1Edges10 = s10.getSequenceFromMask(mask=('[#8 ]', ), )
# region2=a.Surface(side1Edges=side1Edges1+side1Edges2+side1Edges3+side1Edges4+\
#     side1Edges5+side1Edges6+side1Edges7+side1Edges8+side1Edges9+\
#     side1Edges10, name='s_Surf-1')
mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-1', 
    createStepName='Initial', master=regionMaster, slave=regionSlave, 
    sliding=FINITE, thickness=ON, interactionProperty='IntProp-1', 
    adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, 
    clearanceRegion=None)

# session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
#      constraints=OFF, connectors=OFF, engineeringFeatures=OFF)

mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
    numDomains=4, numGPUs=0)


