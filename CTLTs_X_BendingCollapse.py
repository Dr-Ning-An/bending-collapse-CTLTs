# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import os


#############################################################################################################################
# Created by Ning An
# 2019/12
# http://www.anning003.com/create-virtual-nodes
# Created in Abaqus Version 2017

# Function for creating virtual nodes
# mdb: model database
# NameModel: A string with the name of your model
# NameRef1 & NameRef2: A string with the name of your two virtual nodes.
#############################################################################################################################
def VirtualNodes(mdb, NameModel, NameRef1, x, y, z):
    from part import THREE_D, DEFORMABLE_BODY
    #Create reference parts and assemble
    mdb.models[NameModel].Part(dimensionality=THREE_D, name=NameRef1, type=
        DEFORMABLE_BODY)
    mdb.models[NameModel].parts[NameRef1].ReferencePoint(point=(x,y,z))
    mdb.models[NameModel].rootAssembly.Instance(dependent=ON, name=NameRef1, 
        part=mdb.models[NameModel].parts[NameRef1])


    #Create set of reference points
    mdb.models[NameModel].rootAssembly.Set(name=NameRef1, referencePoints=(
        mdb.models[NameModel].rootAssembly.instances[NameRef1].referencePoints[1],))


Mdb()

h_total = 8.0 + 10.6*105.0/180.0*pi
meshSize = h_total/50.0

l = 504.0 # Length
w = 1.0 # Web
L = 3.0 # Lumbus
r = (h_total - 2.0*w - L)*3.0/4.0/pi  # Radius

# Loading
Load_Angle = 0.2



session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

FreqModel = 'FreqModel'
mdb.models.changeKey(fromName='Model-1', toName=FreqModel)

mdb.models[FreqModel].ConstrainedSketch(name='__profile__', sheetSize=400.0)
mdb.models[FreqModel].sketches['__profile__'].Line(
    point1=(0.0, 0.0), 
    point2=(w, 0.0))
mdb.models[FreqModel].sketches['__profile__'].ArcByCenterEnds(
    center=(w, r), 
    direction=COUNTERCLOCKWISE, 
    point1=(w, 0.0), 
    point2=(w + r*sin(60.0*pi/180.0), r*(1-cos(60.0*pi/180.0))))
mdb.models[FreqModel].sketches['__profile__'].ArcByCenterEnds(
    center=(w + r*sin(60.0*pi/180.0)*2.0, 0.0), 
    direction=CLOCKWISE, 
    point1=(w + r*sin(60.0*pi/180.0), r*(1-cos(60.0*pi/180.0))), 
    point2=(w + r*sin(60.0*pi/180.0)*2.0, r))
mdb.models[FreqModel].sketches['__profile__'].Line(
    point1=(w + r*sin(60.0*pi/180.0)*2.0, r), 
    point2=(w + r*sin(60.0*pi/180.0)*2.0 + L/2.0, r))

mdb.models[FreqModel].sketches['__profile__'].ConstructionLine(
	point1=(w + r*sin(60.0*pi/180.0)*2.0 + L/2.0, r), point2=(w + r*sin(60.0*pi/180.0)*2.0 + L/2.0, r/2.0))
mdb.models[FreqModel].sketches['__profile__'].copyMirror(
    mirrorLine=mdb.models[FreqModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0, 15.0), ), 
    objectList=(
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt((w/2.0, 0.0), ), 
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt((w + r*sin(30.0*pi/180.0), r*(1-cos(30.0*pi/180.0))), ), 
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0), r*sin(60.0*pi/180.0)), ), 
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0, r), )))

mdb.models[FreqModel].sketches['__profile__'].ConstructionLine(
    point1=(0.0, 0.0), point2=(10.0, 0.0))
mdb.models[FreqModel].sketches['__profile__'].copyMirror(
    mirrorLine=mdb.models[FreqModel].sketches['__profile__'].geometry.findAt((10.0, 0.0), ), 
    objectList=(
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt((w + r*sin(30.0*pi/180.0), r*(1-cos(30.0*pi/180.0))), ), 
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0), r*sin(60.0*pi/180.0)), ), 
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0, r), ), 
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt(((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(30.0*pi/180.0)), r*(1-cos(30.0*pi/180.0))), ), 
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt(((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0)), r*sin(60.0*pi/180.0)), ), 
    mdb.models[FreqModel].sketches['__profile__'].geometry.findAt(((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0), r), )))


mdb.models[FreqModel].ConstrainedSketch(name='Sketch-1', objectToCopy=
    mdb.models[FreqModel].sketches['__profile__'])

mdb.models[FreqModel].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models[FreqModel].parts['Part-1'].BaseShellExtrude(depth=l, sketch=
    mdb.models[FreqModel].sketches['__profile__'])


mdb.models[FreqModel].parts['Part-1'].Set(faces=
    mdb.models[FreqModel].parts['Part-1'].faces.findAt(((w/2.0, 0.0, 0.0), ), 
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - w/2.0, 0.0, 0.0), ), 
                                                       ), name='Set-Web')
mdb.models[FreqModel].parts['Part-1'].Set(faces=
    mdb.models[FreqModel].parts['Part-1'].faces.findAt(((w + r*sin(30.0*pi/180.0), r*(1-cos(30.0*pi/180.0)), 0.0), ), 
                                                       ((w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0), r*sin(60.0*pi/180.0), 0.0), ), 
                                                       ((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0, r, 0.0), ), 
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(30.0*pi/180.0)), r*(1-cos(30.0*pi/180.0)), 0.0), ),
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0)), r*sin(60.0*pi/180.0), 0.0), ),
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0), r, 0.0), ), 
                                                       ((w + r*sin(30.0*pi/180.0), -r*(1-cos(30.0*pi/180.0)), 0.0), ), 
                                                       ((w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0), -r*sin(60.0*pi/180.0), 0.0), ), 
                                                       ((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0, -r, 0.0), ), 
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(30.0*pi/180.0)), -r*(1-cos(30.0*pi/180.0)), 0.0), ),
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 - r*cos(60.0*pi/180.0)), -r*sin(60.0*pi/180.0), 0.0), ),
                                                       (((w + r*sin(60.0*pi/180.0)*2.0 + L/2.0)*2.0 - (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0/2.0), -r, 0.0), ), 
                                                       ), name='Set-Flanges')

mdb.models[FreqModel].Material(name='Material-1')
mdb.models['FreqModel'].materials['Material-1'].Density(table=((1.2e-09, ), 
    ))
mdb.models[FreqModel].materials['Material-1'].Elastic(
    table=((128.0e3, 6.5e3, 0.35, 7.5e3, 7.5e3, 7.5e3), ), type=LAMINA)

mdb.models[FreqModel].parts['Part-1'].MaterialOrientation(
    additionalRotationType=ROTATION_NONE, axis=AXIS_2, fieldName='', localCsys=
    None, orientationType=GLOBAL, region=
    mdb.models[FreqModel].parts['Part-1'].sets['Set-Flanges'])
mdb.models[FreqModel].parts['Part-1'].MaterialOrientation(
    additionalRotationType=ROTATION_NONE, axis=AXIS_2, fieldName='', localCsys=
    None, orientationType=GLOBAL, region=
    mdb.models[FreqModel].parts['Part-1'].sets['Set-Web'])

t = 0.071
mdb.models[FreqModel].CompositeShellSection(
    idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, 
    layup=(SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1'), 
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'),  
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'), 
           SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1')), 
    name='Section-Flanges', 
    poissonDefinition=DEFAULT, preIntegrate=OFF, symmetric=False, temperature=GRADIENT, thicknessModulus=
    None, thicknessType=UNIFORM, useDensity=OFF)

mdb.models[FreqModel].CompositeShellSection(
    idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, 
    layup=(SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1'), 
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'),  
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'), 
           SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1'),
           SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1'),
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'),
           SectionLayer(thickness=t/4, orientAngle=90.0, material='Material-1'),
           SectionLayer(thickness=t/4, orientAngle=0.0, material='Material-1')), 
    name='Section-Web', 
    poissonDefinition=DEFAULT, preIntegrate=OFF, symmetric=False, temperature=GRADIENT, thicknessModulus=
    None, thicknessType=UNIFORM, useDensity=OFF)


mdb.models[FreqModel].parts['Part-1'].SectionAssignment(
    offset=0.0, 
    offsetField='', 
    offsetType=MIDDLE_SURFACE, 
    region=mdb.models[FreqModel].parts['Part-1'].sets['Set-Web'], 
    sectionName='Section-Web', 
    thicknessAssignment=FROM_SECTION)


mdb.models[FreqModel].parts['Part-1'].SectionAssignment(
    offset=0.0, 
    offsetField='', 
    offsetType=BOTTOM_SURFACE, 
    region=mdb.models[FreqModel].parts['Part-1'].sets['Set-Flanges'], 
    sectionName='Section-Flanges', 
    thicknessAssignment=FROM_SECTION)

mdb.models[FreqModel].rootAssembly.Instance(dependent=OFF, name='Part-1-1', 
    part=mdb.models[FreqModel].parts['Part-1'])

mdb.models[FreqModel].rootAssembly.seedPartInstance(deviationFactor=0.1, 
    minSizeFactor=0.1, regions=(
    mdb.models[FreqModel].rootAssembly.instances['Part-1-1'], ), size=meshSize)
mdb.models[FreqModel].rootAssembly.generateMesh(regions=(
    mdb.models[FreqModel].rootAssembly.instances['Part-1-1'], ))

Bottom_Nodes = mdb.models[FreqModel].rootAssembly.instances['Part-1-1'].nodes.getByBoundingBox(-200.0, -200.0, -0.0001, 200.0, 200.0, 0.0001) 
mdb.models[FreqModel].rootAssembly.Set(name='Bottom_Nodes', nodes=Bottom_Nodes)
Top_Nodes = mdb.models[FreqModel].rootAssembly.instances['Part-1-1'].nodes.getByBoundingBox(-200.0, -200.0, l-0.0001, 200.0, 200.0, l+0.0001) 
mdb.models[FreqModel].rootAssembly.Set(name='Top_Nodes', nodes=Top_Nodes)

# # # #
VirtualNodes(mdb, FreqModel, 'Ref-1', (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0), 0.0, 0.0)
VirtualNodes(mdb, FreqModel, 'Ref-2', (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0), 0.0, l)
VirtualNodes(mdb, FreqModel, 'Ref-0', (w + r*sin(60.0*pi/180.0)*2.0 + L/2.0), 0.0, l/2.0)
mdb.models[FreqModel].Coupling(controlPoint=
    mdb.models[FreqModel].rootAssembly.sets['Ref-1'], couplingType=KINEMATIC, 
    influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-1', 
    surface=mdb.models[FreqModel].rootAssembly.sets['Bottom_Nodes'], u1=ON, u2=ON, 
    u3=ON, ur1=ON, ur2=ON, ur3=ON)
mdb.models[FreqModel].Coupling(controlPoint=
    mdb.models[FreqModel].rootAssembly.sets['Ref-2'], couplingType=KINEMATIC, 
    influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-2', 
    surface=mdb.models[FreqModel].rootAssembly.sets['Top_Nodes'], u1=ON, u2=ON, 
    u3=ON, ur1=ON, ur2=ON, ur3=ON)
mdb.models[FreqModel].Equation(name='Constraint-3', terms=((1.0, 'Ref-2', 
    4), (-1.0, 'Ref-1', 4), (-1.0, 'Ref-0', 4)))
mdb.models[FreqModel].Equation(name='Constraint-4', terms=((1.0, 'Ref-2', 
    3), (1.0, 'Ref-1', 3), (-1.0, 'Ref-0', 3)))


mdb.models[FreqModel].FrequencyStep(name='FreqStep', numEigen=10, previous=
    'Initial')

# # # #
mdb.models[FreqModel].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
    region=mdb.models[FreqModel].rootAssembly.sets['Ref-1'], u1=SET, u2=SET, 
    u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET)
mdb.models[FreqModel].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
    region=mdb.models[FreqModel].rootAssembly.sets['Ref-2'], u1=SET, u2=SET, 
    u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET)


jobName_Freq = 'CTLTs_X_Freq'
mdb.Job(model=FreqModel, name=jobName_Freq)
mdb.saveAs(pathName=jobName_Freq)
mdb.jobs[jobName_Freq].writeInput(consistencyChecking=OFF)
mdb.jobs[jobName_Freq].submit()
mdb.jobs[jobName_Freq].waitForCompletion()


#######################################################
## Extract Eigen Frequencies
#######################################################
jobName = jobName_Freq
stepName = 'FreqStep'
from odbAccess import *
from abaqusConstants import *
import string
import numpy as np
import os
odb = openOdb(path = jobName + '.odb')
freqData = odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['EIGFREQ'].data
np.savetxt(jobName + '.csv', freqData, fmt='%6.12f', delimiter=',')
odb.close()

FundamentalFreq = freqData[1][1]
ExplicitTimePeriod = 30.0*1.0/FundamentalFreq
eta = 1.0E-6


#=========================================================================
# Delete unnecessary steps
#=========================================================================
#---------------------------------------------------------------
ExplicitModel = 'ExplicitModel'
mdb.Model(name=ExplicitModel, objectToCopy=mdb.models[FreqModel])
del mdb.models[ExplicitModel].steps['FreqStep']

print(FundamentalFreq)

mdb.models[ExplicitModel].rootAssembly.setElementType(elemTypes=(ElemType(
    elemCode=S4R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, 
    hourglassControl=ENHANCED), ), 
    regions=(mdb.models[ExplicitModel].rootAssembly.instances['Part-1-1'].faces, ))

mdb.models[ExplicitModel].ContactProperty('IntProp-1')
mdb.models[ExplicitModel].interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=FRICTIONLESS)
mdb.models[ExplicitModel].interactionProperties['IntProp-1'].NormalBehavior(
    allowSeparation=ON, constraintEnforcementMethod=DEFAULT, 
    pressureOverclosure=HARD)
mdb.models[ExplicitModel].ContactExp(createStepName='Initial', name='Int-1')
mdb.models[ExplicitModel].interactions['Int-1'].includedPairs.setValuesInStep(
    stepName='Initial', useAllstar=ON)
mdb.models[ExplicitModel].interactions['Int-1'].contactPropertyAssignments.appendInStep(
    assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Initial')

mdb.models[ExplicitModel].ExplicitDynamicsStep(improvedDtMethod=ON, 
    name='Step-1', previous='Initial', timePeriod=ExplicitTimePeriod)
mdb.models[ExplicitModel].fieldOutputRequests['F-Output-1'].setValues(
    timeInterval=ExplicitTimePeriod/1000)
mdb.models[ExplicitModel].historyOutputRequests['H-Output-1'].setValues(
    timeInterval=ExplicitTimePeriod/1000)
mdb.models[ExplicitModel].fieldOutputRequests['F-Output-1'].setValues(
    variables=('U', 'RF'))
mdb.models['ExplicitModel'].historyOutputRequests['H-Output-1'].setValues(
    variables=('ALLAE', 'ALLIE', 'ALLKE', 'ALLSE', 'ALLVD', 'ALLWK', 'ETOTAL'))
# # # #
mdb.models[ExplicitModel].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
    region=mdb.models[ExplicitModel].rootAssembly.sets['Ref-1'], u1=SET, u2=SET, 
    u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET)
mdb.models[ExplicitModel].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
    region=mdb.models[ExplicitModel].rootAssembly.sets['Ref-2'], u1=SET, u2=SET, 
    u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET)

mdb.models[ExplicitModel].SmoothStepAmplitude(data=((0.0, 0.0), (ExplicitTimePeriod, 
    1.0)), name='Amp-1', timeSpan=STEP)
mdb.models['ExplicitModel'].DisplacementBC(amplitude='Amp-1', createStepName=
    'Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None
    , name='BC-3', region=
    mdb.models['ExplicitModel'].rootAssembly.sets['Ref-0'], u1=SET, u2=SET, 
    u3=SET, ur1=-Load_Angle, ur2=SET, ur3=SET)

mdb.models['ExplicitModel'].rootAssembly.Surface(name='Surf-1', side1Faces=
    mdb.models['ExplicitModel'].rootAssembly.instances['Part-1-1'].faces)
mdb.models['ExplicitModel'].Pressure(amplitude=UNSET, createStepName='Step-1', 
    distributionType=VISCOUS, field='', magnitude=eta*(1.2e-09*128.0e3)**0.5, name='Load-1', region=
    mdb.models['ExplicitModel'].rootAssembly.surfaces['Surf-1'])

jobName_Explicit = 'CTLTs_X_Explicit'
mdb.Job(model=ExplicitModel, name=jobName_Explicit, numCpus=12, numDomains=12, numGPUs=1, explicitPrecision=
    DOUBLE_PLUS_PACK, nodalOutputPrecision=FULL)
mdb.saveAs(pathName=jobName_Explicit)
mdb.jobs[jobName_Explicit].writeInput(consistencyChecking=OFF)
mdb.jobs[jobName_Explicit].submit()
mdb.jobs[jobName_Explicit].waitForCompletion()

##############################################################
## Report Displacement and Reaction Force
## http://www.anning003.com/extract-reaction-force/
##############################################################

jobName = jobName_Explicit
stepName = 'Step-1'
outputSetName = 'REF-0'

from odbAccess import*
from abaqusConstants import*
import string
import numpy as np
import os

odb = openOdb(path = jobName + '.odb')

outfile = open(jobName + '.csv', 'w')
outfile.write('Displacement UR1 [rad]' + ',' + 'Reaction Moment RM1 [Nm]' + '\n')

for fm in range(0, len(odb.steps[stepName].frames)):

  timeFrame = odb.steps[stepName].frames[fm]
  readNode = odb.rootAssembly.nodeSets[outputSetName.upper()]
  Disp = timeFrame.fieldOutputs['UR']
  ReForce = timeFrame.fieldOutputs['RM']
  readNodeDisp = Disp.getSubset(region=readNode)
  readNodeDispValues = readNodeDisp.values
  readNodeRF = ReForce.getSubset(region=readNode)
  readNodeRFValues = readNodeRF.values

  Displacement = np.zeros(len(odb.steps[stepName].frames))
  ReactionForce = np.zeros(len(odb.steps[stepName].frames))
  Displacement[fm] = readNodeDispValues[0].dataDouble[0] # 0-X Direction; 1-Y Direction; 2-Z Direction
  ReactionForce[fm] = readNodeRFValues[0].dataDouble[0]/1000.0

  outfile.write(str(-Displacement[fm]) + ',' + str(-ReactionForce[fm]) + ',' + '\n')

outfile.close()

outfile = open(jobName + '_Energy.csv', 'w')
outfile.write('Time [s]' + ',' + 'EI [mJ]' + ',' + 'EE [mJ]' + ',' + 'EA [mJ]' + ',' + \
    'EV [mJ]' + ',' + 'EKE [mJ]' + ',' + 'EW [mJ]' + ',' + 'ETOT [mJ]' + '\n')

n = len(odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLIE'].data)
timeInc = np.zeros(n)
EI, EE, EA, EV, EKE, EW, ETOT = [np.zeros(n) for i in range(7)]
# EI: Internal Energy   EE: Strain Energy   EA: Artificial Strain Energy
# EV: Viscous Disspated Energy  EKE: Kinetic Energy
# EW: Work Done by External Forces
# ETOT: Energy Balance

for i in range(0, n):

    timeInc[i] = odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLIE'].data[i][0]
    EI[i] = odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLIE'].data[i][1]
    EE[i] = odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLSE'].data[i][1]
    EA[i] = odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLAE'].data[i][1]
    EV[i] = odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLVD'].data[i][1]
    EKE[i] = odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLKE'].data[i][1]
    EW[i] = odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLWK'].data[i][1]
    ETOT[i] = odb.steps[stepName].historyRegions['Assembly ASSEMBLY'].historyOutputs['ETOTAL'].data[i][1]

    outfile.write(str(timeInc[i]) + ',' + str(EI[i]) + ',' + str(EE[i]) + ',' + str(EA[i]) + \
        ',' + str(EV[i]) + ',' + str(EKE[i]) + ',' + str(EW[i]) + ',' + str(ETOT[i]) + ',' + '\n')

outfile.close()


odb.close()