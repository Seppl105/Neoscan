"""
modelAExample.py

A simple example: Creating a part.
"""

from abaqus import *
from abaqusConstants import *
import sketch
import part

myModel = mdb.Model(name='Model A')

mySketch = myModel.Sketch(name='Sketch A', sheetSize=200.0)

xyCoordsInner = ((-5 , 20), (5, 20), (15, 0),
    (-15, 0), (-5, 20))

xyCoordsOuter = ((-10, 30), (10, 30), (40, -30),
    (30, -30), (20, -10), (-20, -10),
    (-30, -30), (-40, -30), (-10, 30))

for i in range(len(xyCoordsInner)-1):
    mySketch.Line(point1=xyCoordsInner[i],
        point2=xyCoordsInner[i+1])

for i in range(len(xyCoordsOuter)-1):
    mySketch.Line(point1=xyCoordsOuter[i],
        point2=xyCoordsOuter[i+1])

myPart = myModel.Part(name='Part A', dimensionality=THREE_D,
    type=DEFORMABLE_BODY)

myPart.BaseSolidExtrude(sketch=mySketch, depth=20.0)

myViewport = session.Viewport(name='Viewport for Model A',
    origin=(10, 10), width=150, height=100)

myViewport.setValues(displayedObject=myPart)

myViewport.partDisplay.setValues(renderStyle=SHADED)
