import math
from odbAccess import *
from abaqusConstants import *
InputName='test1'
print 'ODB = test1'
outputfilename=InputName+'-Modulii.txt'
out = open(outputfilename,'w')
odb = openOdb(InputName+'.odb')
numberFrame = len(odb.steps['STEP-1'].frames)
nset=odb.rootAssembly.instances['PART-1-1'].nodeSets['RIGHTSURFACE']
for ns in range (1,numberFrame):
	frame = odb.steps['STEP-1'].frames[ns]
	force=frame.fieldOutputs['RF']
	nsetForce = force.getSubset(region=nset)
	sumforceReal=0
	sumforceImag=0
	for v in nsetForce.values:
		sumforceReal+=v.data[0]/(0.005*4*4/4)
		sumforceImag+=v.conjugateData[0]/(0.005*4*4/4)
	out.write(str(frame.frameValue*2*math.pi)+'\t'+str(sumforceReal)+'\t'+str(sumforceImag)+'\n')
out.close()
odb.close()
