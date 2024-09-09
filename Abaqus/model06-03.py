# coding=utf-8
from abaqus import *
from abaqusConstants import *
import regionToolset

#####################C:\Users\SGnir\OneDrive\TU\Bachelorsemester06\Kernspintomographie (MRT)\GitHub-Neoscan\Neoscan\Abaqus\model06.py
    

# the positive x axis coincides with the r-direction in the modeled plane; the revolution wil be arounnd the y-axis

roundToDigit = 9 # if a value is rounded it will be to the "roundToDigit" digits

# windingDivisor = 60       # Es wird für 600/windingDivisor Windungen gerechnet
# totalWindings = int(600/windingDivisor)
totalWindings = 600 #600
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
yCoorEnd = 0.004            # [m] y Coordinate am Ende des ersten Pancakes
yCoorMid = (yCoorEnd + yCoorStart)/2

j = 114.2 * 10**6                        # [A/m^2] Stromdichte
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

# Create the first ring
part = createRing(name=parts[0] + "1", pointOne=(r_i, yCoorStart), pointTwo=(r_i + materialWidths[0], yCoorEnd))
# assign section
faceOfPart = part.faces.findAt(((r_i + materialWidths[0]/2, yCoorMid, 0.0), ))  # the argument must be tuple of the shape ((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),)
regionOfSectionName= regionToolset.Region(faces=faceOfPart)
part.SectionAssignment(region=regionOfSectionName, sectionName=sectionNames[0])

# add the first ring to the assembly
assembly = mdb.models['Model-1'].rootAssembly
assembly.Instance(name="AssemblyInstance-" + parts[0] + "1", part=part, dependent=ON)

# sectionPoints = [[] for element in sectionNames] # a list containing one list for each section with each list containing one point within each corresponding face
#                                                  # ("sectionPoints" is later used to assign the sections)                                             

# Create all remaining parts, add them to the assembly and update sectionPoints
currentR = r_i    # outer radius of the partioning (rises which each added layer)
for winding in range(totalWindings):
    for i,width in enumerate(materialWidths):
        if winding != 0 or i != 0:
            # Create ring
            part = createRing(name=parts[i] + str(winding + 1), pointOne=(currentR, yCoorStart), pointTwo=(currentR + width, yCoorEnd))
            # Add ring to the assembly
            instanceName = "AssemblyInstance" + parts[i] + str(winding + 1)
            assembly.Instance(name=instanceName, part=part, dependent=ON)
            
            # assign section
            faceOfPart = part.faces.findAt(((currentR + width/2, yCoorMid, 0.0), ))  # the argument must be tuple of the shape ((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),)
            regionOfSectionName= regionToolset.Region(faces=faceOfPart)
            part.SectionAssignment(region=regionOfSectionName, sectionName=sectionNames[i])
            # old:
            # # assign a point within each instance to the corresponding list (used to assign sections) 
            # sectionPoints[i].append((currentR + width/2, yCoorMid, 0.0))
            currentR = round(currentR + width, roundToDigit)
        else:
            currentR = currentR + width
            
            
            

print("currentR equal to assigned or calculated r_a: ", currentR == r_a, " currentR: ", currentR, " r_a: ", r_a)
assembly.regenerate()

# #assign the faces to the corresponding sections according to sectionPoints
# for i,sectionName in enumerate(sectionNames):
#     tuplePoints = tuple([(point,) for point in sectionPoints[i]]) # creating a tuple of tuple ==> tuplPoints => (((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),))
#     facesOfSectionName = part.faces.findAt(*tuplePoints)  # ==> the argument must be tuple of the shape ((0, 0, 0),), ((1, 1, 1),), ((2, 2, 2),) ==> hence we must remove the outer tuple ==> hence, we passed the argument with '*' sign which takes care of this
    
#     #faceOfNewLayerSequence = partPancake.faces.sequenceFromLabels(labels=[faceOfNewLayer.index])
#     regionOfSectionName= regionToolset.Region(faces=facesOfSectionName)
#     partPancake.SectionAssignment(region=regionOfSectionName, sectionName=sectionName)