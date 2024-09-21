# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

def CoilDemo():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(40.0, 5.0))
    mdb.models['Model-1'].ConstrainedSketch(name='Sketch-Coil', objectToCopy=s)
    mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
        toName='Sketch-Coil')
    s.unsetPrimaryObject()
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.sketchOptions.setValues(viewStyle=AXISYM)
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.FixedConstraint(entity=g[2])
    s1.rectangle(point1=(5.0, 0.0), point2=(55.0, 5.0))
    mdb.models['Model-1'].ConstrainedSketch(name='Sketch-2', objectToCopy=s1)
    p = mdb.models['Model-1'].Part(name='Part-Coil', dimensionality=AXISYMMETRIC, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-Coil']
    p.BaseShell(sketch=s1)
    s1.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-Coil']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']


def PancakeDemo():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.sketchOptions.setValues(viewStyle=AXISYM)
    s.setPrimaryObject(option=STANDALONE)
    s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s.FixedConstraint(entity=g[2])
    s.rectangle(point1=(10.0, 5.0), point2=(50.0, 0.0))
    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=AXISYMMETRIC, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-1']
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    
def SectionAssign():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    region = p.Set(faces=faces, name='Set-Assign')
    p = mdb.models['Model-1'].parts['Part-1']
    p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)


def SectionAss():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    region = p.Set(faces=faces, name='Set-4')
    p = mdb.models['Model-1'].parts['Part-1']
    p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)


def Assembly():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    a1 = mdb.models['Model-1'].rootAssembly
    p = mdb.models['Model-1'].parts['Part-1']
    a1.Instance(name='Part-1-1', part=p, dependent=OFF)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=50.8187, 
        farPlane=87.8031, width=57.7306, height=41.81, viewOffsetX=19.5512, 
        viewOffsetY=1.50981)
    a1 = mdb.models['Model-1'].rootAssembly
    a1.LinearInstancePattern(instanceList=('Part-1-1', ), direction1=(1.0, 0.0, 
        0.0), direction2=(0.0, 1.0, 0.0), number1=4, number2=1, spacing1=30.0, 
        spacing2=20.0)


def ContactBC():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON, 
        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    mdb.models['Model-1'].ContactProperty('IntProp-1')
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((0.9, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    mdb.models['Model-1'].ContactExp(name='Int-1', createStepName='Initial')
    mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
        stepName='Initial', useAllstar=ON)
    mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
        stepName='Initial', assignments=((GLOBAL, SELF, 'IntProp-1'), ))
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, interactions=OFF, constraints=OFF, 
        engineeringFeatures=OFF)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.862925, 
        farPlane=0.864381, width=0.00419459, height=0.00172063, 
        viewOffsetX=0.215163, viewOffsetY=-0.00201686)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstance-Con1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#1 ]', ), )
    v1 = a.instances['AssemblyInstance-Con1'].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#3 ]', ), )
    e2 = a.instances['AssemblyInstanceCow1'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#1 ]', ), )
    v2 = a.instances['AssemblyInstanceCow1'].vertices
    verts2 = v2.getSequenceFromMask(mask=('[#3 ]', ), )
    e3 = a.instances['AssemblyInstanceIns1'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#1 ]', ), )
    v3 = a.instances['AssemblyInstanceIns1'].vertices
    verts3 = v3.getSequenceFromMask(mask=('[#3 ]', ), )
    e4 = a.instances['AssemblyInstanceCon2'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#1 ]', ), )
    v4 = a.instances['AssemblyInstanceCon2'].vertices
    verts4 = v4.getSequenceFromMask(mask=('[#3 ]', ), )
    e5 = a.instances['AssemblyInstanceCow2'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#1 ]', ), )
    v5 = a.instances['AssemblyInstanceCow2'].vertices
    verts5 = v5.getSequenceFromMask(mask=('[#3 ]', ), )
    e6 = a.instances['AssemblyInstanceIns2'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#1 ]', ), )
    v6 = a.instances['AssemblyInstanceIns2'].vertices
    verts6 = v6.getSequenceFromMask(mask=('[#3 ]', ), )
    e7 = a.instances['AssemblyInstanceCon3'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#1 ]', ), )
    v7 = a.instances['AssemblyInstanceCon3'].vertices
    verts7 = v7.getSequenceFromMask(mask=('[#3 ]', ), )
    e8 = a.instances['AssemblyInstanceCow3'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#1 ]', ), )
    v8 = a.instances['AssemblyInstanceCow3'].vertices
    verts8 = v8.getSequenceFromMask(mask=('[#3 ]', ), )
    e9 = a.instances['AssemblyInstanceIns3'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#1 ]', ), )
    v9 = a.instances['AssemblyInstanceIns3'].vertices
    verts9 = v9.getSequenceFromMask(mask=('[#3 ]', ), )
    e10 = a.instances['AssemblyInstanceCon4'].edges
    edges10 = e10.getSequenceFromMask(mask=('[#1 ]', ), )
    v10 = a.instances['AssemblyInstanceCon4'].vertices
    verts10 = v10.getSequenceFromMask(mask=('[#3 ]', ), )
    e11 = a.instances['AssemblyInstanceCow4'].edges
    edges11 = e11.getSequenceFromMask(mask=('[#1 ]', ), )
    v11 = a.instances['AssemblyInstanceCow4'].vertices
    verts11 = v11.getSequenceFromMask(mask=('[#3 ]', ), )
    e12 = a.instances['AssemblyInstanceIns4'].edges
    edges12 = e12.getSequenceFromMask(mask=('[#1 ]', ), )
    v12 = a.instances['AssemblyInstanceIns4'].vertices
    verts12 = v12.getSequenceFromMask(mask=('[#3 ]', ), )
    e13 = a.instances['AssemblyInstanceCon5'].edges
    edges13 = e13.getSequenceFromMask(mask=('[#1 ]', ), )
    v13 = a.instances['AssemblyInstanceCon5'].vertices
    verts13 = v13.getSequenceFromMask(mask=('[#3 ]', ), )
    e14 = a.instances['AssemblyInstanceCow5'].edges
    edges14 = e14.getSequenceFromMask(mask=('[#1 ]', ), )
    v14 = a.instances['AssemblyInstanceCow5'].vertices
    verts14 = v14.getSequenceFromMask(mask=('[#3 ]', ), )
    e15 = a.instances['AssemblyInstanceIns5'].edges
    edges15 = e15.getSequenceFromMask(mask=('[#1 ]', ), )
    v15 = a.instances['AssemblyInstanceIns5'].vertices
    verts15 = v15.getSequenceFromMask(mask=('[#3 ]', ), )
    region = a.Set(vertices=verts1+verts2+verts3+verts4+verts5+verts6+verts7+\
        verts8+verts9+verts10+verts11+verts12+verts13+verts14+verts15, 
        edges=edges1+edges2+edges3+edges4+edges5+edges6+edges7+edges8+edges9+\
        edges10+edges11+edges12+edges13+edges14+edges15, name='Set-BC')
    mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Initial', 
        region=region, u1=UNSET, u2=SET, ur3=UNSET, amplitude=UNSET, 
        distributionType=UNIFORM, fieldName='', localCsys=None)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
        predefinedFields=OFF, connectors=OFF, adaptiveMeshConstraints=ON)
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')


def StepContactAgain():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON, 
        adaptiveMeshConstraints=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863209, 
        farPlane=0.864097, width=0.00289467, height=0.0011874, 
        viewOffsetX=0.215196, viewOffsetY=-0.00197148)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
        leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)


def Interaction():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    leaf = dgm.LeafFromGeometry(edgeSeq=edges1+edges2+edges3+edges4+edges5+edges6+\
        edges7+edges8+edges9)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    leaf = dgm.LeafFromGeometry(edgeSeq=edges1+edges2+edges3+edges4+edges5+edges6+\
        edges7+edges8+edges9)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    leaf = dgm.LeafFromGeometry(edgeSeq=edges1+edges2+edges3+edges4+edges5+edges6+\
        edges7+edges8+edges9)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863273, 
        farPlane=0.864034, width=0.00219217, height=0.000899233, 
        viewOffsetX=0.214648, viewOffsetY=-0.00177058)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    leaf = dgm.LeafFromGeometry(edgeSeq=edges1+edges2+edges3+edges4+edges5+edges6+\
        edges7+edges8+edges9)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    leaf = dgm.LeafFromGeometry(edgeSeq=edges1+edges2+edges3+edges4+edges5+edges6+\
        edges7+edges8+edges9)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863374, 
        farPlane=0.863933, width=0.00182099, height=0.000746976, 
        viewOffsetX=0.214968, viewOffsetY=-0.00175037)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863249, 
        farPlane=0.864058, width=0.00233203, height=0.000956604, 
        viewOffsetX=0.214949, viewOffsetY=-0.00170251)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    leaf = dgm.LeafFromGeometry(edgeSeq=edges1+edges2+edges3+edges4+edges5+edges6+\
        edges7+edges8+edges9)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
        leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    leaf = dgm.LeafFromGeometry(edgeSeq=edges1+edges2+edges3+edges4+edges5+edges6+\
        edges7+edges8+edges9)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    leaf = dgm.LeafFromGeometry(edgeSeq=edges1+edges2+edges3+edges4+edges5+edges6+\
        edges7+edges8+edges9)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    leaf = dgm.LeafFromGeometry(edgeSeq=edges1+edges2+edges3+edges4+edges5+edges6+\
        edges7+edges8+edges9)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
        leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow2'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon3'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow3'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    e6 = a.instances['AssemblyInstanceCon4'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#8 ]', ), )
    e7 = a.instances['AssemblyInstanceCow4'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#2 ]', ), )
    e8 = a.instances['AssemblyInstanceCon5'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#8 ]', ), )
    e9 = a.instances['AssemblyInstanceCow5'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#2 ]', ), )
    a.Set(edges=edges1+edges2+edges3+edges4+edges5+edges6+edges7+edges8+edges9, 
        name='Set-MasterSurfaces')
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863296, 
        farPlane=0.864011, width=0.00206069, height=0.000845301, 
        viewOffsetX=0.214891, viewOffsetY=-0.00172383)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    f1 = a.instances['AssemblyInstanceCow2'].faces
    faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
    f2 = a.instances['AssemblyInstanceCow1'].faces
    faces2 = f2.getSequenceFromMask(mask=('[#1 ]', ), )
    f3 = a.instances['AssemblyInstanceCow3'].faces
    faces3 = f3.getSequenceFromMask(mask=('[#1 ]', ), )
    f4 = a.instances['AssemblyInstanceCow4'].faces
    faces4 = f4.getSequenceFromMask(mask=('[#1 ]', ), )
    f5 = a.instances['AssemblyInstanceCow5'].faces
    faces5 = f5.getSequenceFromMask(mask=('[#1 ]', ), )
    leaf = dgm.LeafFromGeometry(faceSeq=faces1+faces2+faces3+faces4+faces5)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#8 ]', ), )
    e2 = a.instances['AssemblyInstanceCow2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceCow3'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#8 ]', ), )
    e4 = a.instances['AssemblyInstanceCow4'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceCow5'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#8 ]', ), )
    a.Set(edges=edges1+edges2+edges3+edges4+edges5, 
        name='Set-MasterSurfacesConToCow')
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstance-Con1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
    e2 = a.instances['AssemblyInstanceCon2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#2 ]', ), )
    e3 = a.instances['AssemblyInstanceCon3'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#2 ]', ), )
    e4 = a.instances['AssemblyInstanceCon4'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#2 ]', ), )
    e5 = a.instances['AssemblyInstanceCon5'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#2 ]', ), )
    a.Set(edges=edges1+edges2+edges3+edges4+edges5, name='Set-SlaveCon')
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
        leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    f1 = a.instances['AssemblyInstanceCow1'].faces
    faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
    f2 = a.instances['AssemblyInstanceCow2'].faces
    faces2 = f2.getSequenceFromMask(mask=('[#1 ]', ), )
    f3 = a.instances['AssemblyInstanceCow3'].faces
    faces3 = f3.getSequenceFromMask(mask=('[#1 ]', ), )
    f4 = a.instances['AssemblyInstanceCow4'].faces
    faces4 = f4.getSequenceFromMask(mask=('[#1 ]', ), )
    f5 = a.instances['AssemblyInstanceCow5'].faces
    faces5 = f5.getSequenceFromMask(mask=('[#1 ]', ), )
    leaf = dgm.LeafFromGeometry(faceSeq=faces1+faces2+faces3+faces4+faces5)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceIns1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#8 ]', ), )
    e2 = a.instances['AssemblyInstanceIns2'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#8 ]', ), )
    e3 = a.instances['AssemblyInstanceIns3'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#8 ]', ), )
    e4 = a.instances['AssemblyInstanceIns4'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#8 ]', ), )
    e5 = a.instances['AssemblyInstanceIns5'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#8 ]', ), )
    a.Set(edges=edges1+edges2+edges3+edges4+edges5, name='Set-SlaveInsToCow')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['AssemblyInstanceCow1'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#a ]', ), )
    s2 = a.instances['AssemblyInstanceCow2'].edges
    side1Edges2 = s2.getSequenceFromMask(mask=('[#a ]', ), )
    s3 = a.instances['AssemblyInstanceCow3'].edges
    side1Edges3 = s3.getSequenceFromMask(mask=('[#a ]', ), )
    s4 = a.instances['AssemblyInstanceCow4'].edges
    side1Edges4 = s4.getSequenceFromMask(mask=('[#a ]', ), )
    s5 = a.instances['AssemblyInstanceCow5'].edges
    side1Edges5 = s5.getSequenceFromMask(mask=('[#a ]', ), )
    region1=regionToolset.Region(side1Edges=side1Edges1+side1Edges2+side1Edges3+\
        side1Edges4+side1Edges5)
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['AssemblyInstance-Con1'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#2 ]', ), )
    s2 = a.instances['AssemblyInstanceIns1'].edges
    side1Edges2 = s2.getSequenceFromMask(mask=('[#8 ]', ), )
    s3 = a.instances['AssemblyInstanceCon2'].edges
    side1Edges3 = s3.getSequenceFromMask(mask=('[#2 ]', ), )
    s4 = a.instances['AssemblyInstanceIns2'].edges
    side1Edges4 = s4.getSequenceFromMask(mask=('[#8 ]', ), )
    s5 = a.instances['AssemblyInstanceCon3'].edges
    side1Edges5 = s5.getSequenceFromMask(mask=('[#2 ]', ), )
    s6 = a.instances['AssemblyInstanceIns3'].edges
    side1Edges6 = s6.getSequenceFromMask(mask=('[#8 ]', ), )
    s7 = a.instances['AssemblyInstanceCon4'].edges
    side1Edges7 = s7.getSequenceFromMask(mask=('[#2 ]', ), )
    s8 = a.instances['AssemblyInstanceIns4'].edges
    side1Edges8 = s8.getSequenceFromMask(mask=('[#8 ]', ), )
    s9 = a.instances['AssemblyInstanceCon5'].edges
    side1Edges9 = s9.getSequenceFromMask(mask=('[#2 ]', ), )
    s10 = a.instances['AssemblyInstanceIns5'].edges
    side1Edges10 = s10.getSequenceFromMask(mask=('[#8 ]', ), )
    region2=regionToolset.Region(side1Edges=side1Edges1+side1Edges2+side1Edges3+\
        side1Edges4+side1Edges5+side1Edges6+side1Edges7+side1Edges8+\
        side1Edges9+side1Edges10)
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-2', 
        createStepName='Initial', master=region1, slave=region2, 
        sliding=FINITE, thickness=ON, interactionProperty='IntProp-1', 
        adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, 
        clearanceRegion=None)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
        leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    a = mdb.models['Model-1'].rootAssembly
    f1 = a.instances['AssemblyInstanceCow1'].faces
    faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
    f2 = a.instances['AssemblyInstanceCow2'].faces
    faces2 = f2.getSequenceFromMask(mask=('[#1 ]', ), )
    f3 = a.instances['AssemblyInstanceCow3'].faces
    faces3 = f3.getSequenceFromMask(mask=('[#1 ]', ), )
    f4 = a.instances['AssemblyInstanceCow4'].faces
    faces4 = f4.getSequenceFromMask(mask=('[#1 ]', ), )
    f5 = a.instances['AssemblyInstanceCow5'].faces
    faces5 = f5.getSequenceFromMask(mask=('[#1 ]', ), )
    f6 = a.instances['AssemblyInstance-Con1'].faces
    faces6 = f6.getSequenceFromMask(mask=('[#1 ]', ), )
    f7 = a.instances['AssemblyInstanceCon2'].faces
    faces7 = f7.getSequenceFromMask(mask=('[#1 ]', ), )
    f8 = a.instances['AssemblyInstanceCon3'].faces
    faces8 = f8.getSequenceFromMask(mask=('[#1 ]', ), )
    f9 = a.instances['AssemblyInstanceCon4'].faces
    faces9 = f9.getSequenceFromMask(mask=('[#1 ]', ), )
    f10 = a.instances['AssemblyInstanceCon5'].faces
    faces10 = f10.getSequenceFromMask(mask=('[#1 ]', ), )
    leaf = dgm.LeafFromGeometry(faceSeq=faces1+faces2+faces3+faces4+faces5+faces6+\
        faces7+faces8+faces9+faces10)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['AssemblyInstanceCon2'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#8 ]', ), )
    s2 = a.instances['AssemblyInstanceCon3'].edges
    side1Edges2 = s2.getSequenceFromMask(mask=('[#8 ]', ), )
    s3 = a.instances['AssemblyInstanceCon4'].edges
    side1Edges3 = s3.getSequenceFromMask(mask=('[#8 ]', ), )
    s4 = a.instances['AssemblyInstanceCon5'].edges
    side1Edges4 = s4.getSequenceFromMask(mask=('[#8 ]', ), )
    region1=regionToolset.Region(side1Edges=side1Edges1+side1Edges2+side1Edges3+\
        side1Edges4)
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['AssemblyInstanceIns1'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#2 ]', ), )
    s2 = a.instances['AssemblyInstanceIns2'].edges
    side1Edges2 = s2.getSequenceFromMask(mask=('[#2 ]', ), )
    s3 = a.instances['AssemblyInstanceIns3'].edges
    side1Edges3 = s3.getSequenceFromMask(mask=('[#2 ]', ), )
    s4 = a.instances['AssemblyInstanceIns4'].edges
    side1Edges4 = s4.getSequenceFromMask(mask=('[#2 ]', ), )
    region2=regionToolset.Region(side1Edges=side1Edges1+side1Edges2+side1Edges3+\
        side1Edges4)
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-3', 
        createStepName='Initial', master=region1, slave=region2, 
        sliding=FINITE, thickness=ON, interactionProperty='IntProp-1', 
        adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, 
        clearanceRegion=None)


def CreateSetForBC():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstance-Con1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#1 ]', ), )
    e2 = a.instances['AssemblyInstanceCow1'].edges
    edges2 = e2.getSequenceFromMask(mask=('[#1 ]', ), )
    e3 = a.instances['AssemblyInstanceIns1'].edges
    edges3 = e3.getSequenceFromMask(mask=('[#1 ]', ), )
    e4 = a.instances['AssemblyInstanceCon2'].edges
    edges4 = e4.getSequenceFromMask(mask=('[#1 ]', ), )
    e5 = a.instances['AssemblyInstanceCow2'].edges
    edges5 = e5.getSequenceFromMask(mask=('[#1 ]', ), )
    e6 = a.instances['AssemblyInstanceIns2'].edges
    edges6 = e6.getSequenceFromMask(mask=('[#1 ]', ), )
    e7 = a.instances['AssemblyInstanceCon3'].edges
    edges7 = e7.getSequenceFromMask(mask=('[#1 ]', ), )
    e8 = a.instances['AssemblyInstanceCow3'].edges
    edges8 = e8.getSequenceFromMask(mask=('[#1 ]', ), )
    e9 = a.instances['AssemblyInstanceIns3'].edges
    edges9 = e9.getSequenceFromMask(mask=('[#1 ]', ), )
    e10 = a.instances['AssemblyInstanceCon4'].edges
    edges10 = e10.getSequenceFromMask(mask=('[#1 ]', ), )
    e11 = a.instances['AssemblyInstanceCow4'].edges
    edges11 = e11.getSequenceFromMask(mask=('[#1 ]', ), )
    e12 = a.instances['AssemblyInstanceIns4'].edges
    edges12 = e12.getSequenceFromMask(mask=('[#1 ]', ), )
    e13 = a.instances['AssemblyInstanceCon5'].edges
    edges13 = e13.getSequenceFromMask(mask=('[#1 ]', ), )
    e14 = a.instances['AssemblyInstanceCow5'].edges
    edges14 = e14.getSequenceFromMask(mask=('[#1 ]', ), )
    e15 = a.instances['AssemblyInstanceIns5'].edges
    edges15 = e15.getSequenceFromMask(mask=('[#1 ]', ), )
    a.Set(edges=edges1+edges2+edges3+edges4+edges5+edges6+edges7+edges8+edges9+\
        edges10+edges11+edges12+edges13+edges14+edges15, name='Set-BC2')
def Interaction3():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    mdb.models['Model-1'].ContactProperty('IntProp-1')
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((0.9, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    mdb.models['Model-1'].interactionProperties['IntProp-1'].tangentialBehavior.setValues(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((0.9, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.either(leaf=leaf)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863083, 
        farPlane=0.864223, width=0.00411057, height=0.00168617, 
        viewOffsetX=0.214927, viewOffsetY=-0.00120241)
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['AssemblyInstanceCow1'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#a ]', ), )
    s2 = a.instances['AssemblyInstanceCow2'].edges
    side1Edges2 = s2.getSequenceFromMask(mask=('[#a ]', ), )
    s3 = a.instances['AssemblyInstanceCow3'].edges
    side1Edges3 = s3.getSequenceFromMask(mask=('[#a ]', ), )
    s4 = a.instances['AssemblyInstanceCow4'].edges
    side1Edges4 = s4.getSequenceFromMask(mask=('[#a ]', ), )
    s5 = a.instances['AssemblyInstanceCow5'].edges
    side1Edges5 = s5.getSequenceFromMask(mask=('[#a ]', ), )
    region1=a.Surface(side1Edges=side1Edges1+side1Edges2+side1Edges3+side1Edges4+\
        side1Edges5, name='m_Surf-1')
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['AssemblyInstance-Con1'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#2 ]', ), )
    s2 = a.instances['AssemblyInstanceIns1'].edges
    side1Edges2 = s2.getSequenceFromMask(mask=('[#8 ]', ), )
    s3 = a.instances['AssemblyInstanceCon2'].edges
    side1Edges3 = s3.getSequenceFromMask(mask=('[#2 ]', ), )
    s4 = a.instances['AssemblyInstanceIns2'].edges
    side1Edges4 = s4.getSequenceFromMask(mask=('[#8 ]', ), )
    s5 = a.instances['AssemblyInstanceCon3'].edges
    side1Edges5 = s5.getSequenceFromMask(mask=('[#2 ]', ), )
    s6 = a.instances['AssemblyInstanceIns3'].edges
    side1Edges6 = s6.getSequenceFromMask(mask=('[#8 ]', ), )
    s7 = a.instances['AssemblyInstanceCon4'].edges
    side1Edges7 = s7.getSequenceFromMask(mask=('[#2 ]', ), )
    s8 = a.instances['AssemblyInstanceIns4'].edges
    side1Edges8 = s8.getSequenceFromMask(mask=('[#8 ]', ), )
    s9 = a.instances['AssemblyInstanceCon5'].edges
    side1Edges9 = s9.getSequenceFromMask(mask=('[#2 ]', ), )
    s10 = a.instances['AssemblyInstanceIns5'].edges
    side1Edges10 = s10.getSequenceFromMask(mask=('[#8 ]', ), )
    region2=a.Surface(side1Edges=side1Edges1+side1Edges2+side1Edges3+side1Edges4+\
        side1Edges5+side1Edges6+side1Edges7+side1Edges8+side1Edges9+\
        side1Edges10, name='s_Surf-1')
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-1', 
        createStepName='Initial', master=region1, slave=region2, 
        sliding=FINITE, thickness=ON, interactionProperty='IntProp-1', 
        adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, 
        clearanceRegion=None)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
        constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
    mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
        numDomains=4, numGPUs=0)


def NoteSetMeshSubmit():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceCow1'].edges
    e2 = a.instances['AssemblyInstanceCow2'].edges
    e3 = a.instances['AssemblyInstanceCow3'].edges
    e4 = a.instances['AssemblyInstanceCow4'].edges
    e5 = a.instances['AssemblyInstanceCow5'].edges
    pickedEdges = e1.getSequenceFromMask(mask=('[#1 ]', ), )+\
        e2.getSequenceFromMask(mask=('[#1 ]', ), )+e3.getSequenceFromMask(
        mask=('[#1 ]', ), )+e4.getSequenceFromMask(mask=('[#1 ]', ), )+\
        e5.getSequenceFromMask(mask=('[#1 ]', ), )
    a.seedEdgeByNumber(edges=pickedEdges, number=3, constraint=FINER)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863513, 
        farPlane=0.863794, width=0.000808477, height=0.000332556, 
        viewOffsetX=0.214646, viewOffsetY=-0.00194995)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['AssemblyInstanceIns1'].edges
    e2 = a.instances['AssemblyInstanceIns2'].edges
    pickedEdges = e1.getSequenceFromMask(mask=('[#1 ]', ), )+\
        e2.getSequenceFromMask(mask=('[#1 ]', ), )
    a.seedEdgeBySize(edges=pickedEdges, size=3.33e-06, deviationFactor=0.1, 
        minSizeFactor=0.1, constraint=FINER)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863154, 
        farPlane=0.864152, width=0.00324404, height=0.00133439, 
        viewOffsetX=0.215096, viewOffsetY=-0.00211053)
    a = mdb.models['Model-1'].rootAssembly
    partInstances =(a.instances['AssemblyInstance-Con1'], 
        a.instances['AssemblyInstanceCow1'], 
        a.instances['AssemblyInstanceIns1'], 
        a.instances['AssemblyInstanceCon2'], 
        a.instances['AssemblyInstanceCow2'], 
        a.instances['AssemblyInstanceIns2'], 
        a.instances['AssemblyInstanceCon3'], 
        a.instances['AssemblyInstanceCow3'], 
        a.instances['AssemblyInstanceIns3'], 
        a.instances['AssemblyInstanceCon4'], 
        a.instances['AssemblyInstanceCow4'], 
        a.instances['AssemblyInstanceIns4'], 
        a.instances['AssemblyInstanceCon5'], 
        a.instances['AssemblyInstanceCow5'], 
        a.instances['AssemblyInstanceIns5'], )
    a.generateMesh(regions=partInstances)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.863082, 
        farPlane=0.864225, width=0.00371632, height=0.00152866, 
        viewOffsetX=0.215131, viewOffsetY=-0.00209442)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
    session.mdbData.summary()
    o3 = session.openOdb(
        name='C:/Users/SGnir/OneDrive/TU/Bachelorsemester06/Kernspintomographie (MRT)/GitHub-Neoscan/Neoscan/Abaqus/Job-1.odb')
    session.viewports['Viewport: 1'].setValues(displayedObject=o3)
    session.viewports['Viewport: 1'].makeCurrent()
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_DEF, ))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00946832, 
        farPlane=0.0153446, width=0.00431538, height=0.00185307, 
        viewOffsetX=0.000638195, viewOffsetY=-0.000133127)
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(
        COMPONENT, 'S11'), )
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00958645, 
        farPlane=0.0152265, width=0.00362901, height=0.00155833, 
        viewOffsetX=0.00070763, viewOffsetY=-0.000111072)
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(
        COMPONENT, 'S22'), )


def MeshPart():
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    pass


