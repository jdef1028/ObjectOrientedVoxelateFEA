__author__ = "Xiaolin Li, Northwestern University"
__copyright__ = "Advanced Materials Lab, Integrated Design Automation Lab"
__credits__ = ["Xiaolin Li", "L. Catherine Brinson", "Wei Chen"]
__version__ = "v0.0.1"
__maintainer__ = "Xiaolin Li"
__email__ = "xiaolinli2018@u.northwestern.edu"
__status__ = "Developing"

import scipy.io
import math
import os
import shutil
import sys
import numpy


class allClusterBase():
	""" This class defines the common methods and attributes that EVERY 3D microstructure could use"""

	def __isolateClusters__(self, matfile):
		""" __isolateClusters__ is a private method that could be only used for the microstructures with multiple clusters.
			This function calls the matlab code to generate files contains individual clusters
			@para: matfile --> the file name of the 3D microstructure .mat file
			@return: None. But generate a folder with multiple .mat files
			Developed on Jul. 14, 2015.
			"""
			
		connectivity = 18 # select from [6, 18, 26] --> the number of neighboring pixels that are taken as connected
		command = 'matlab -r "'+"cluster_iso('"+ str(matfile) + "',"+ str(connectivity) + ')"' +' -nojvm -nosplash -nodesktop'
		os.system(command)
		print "==========================================================\n"
		print "==Clusters have been isolated into different .mat files!==\n"
		print "==========================================================\n"

	def loadStructure(self, matfile, var_name):
		""" loadStructure() is a public method that could be used for every 3D microstructure provided.
			This function could load the microstructure in .mat file into an 3D array in python to be processed.
			@para: matfile --> the file name of the microstructure
				   var_name --> the name of the variable in the .mat file that contains the 3D microstructure 0-1 matrix.
			@return: None
			Developed on Jul. 14, 2015.
		"""
		Data = scipy.io.loadmat(matfile, matlab_compatible = True) # import data from .mat file
		self.Image = Data[var_name]
		self.Z = len(self.Image) # Z dimension
		self.Y = len(self.Image[0]) # Y dimension
		self.X = len(self.Image[0][0]) # X dimension
		self.VF = self.getVF()
		print "==========================================================\n"
		print "==            Microstructure has been loaded!           ==\n"
		print "==             " + matfile
		print "==========================================================\n"

	def getVF(self):
		""" getVF() is a public method that could calculate the volume fraction from the loaded 3D microstructure 
			@para: None
			@return: Volume fraction of the microstructure
		"""
		count = 0
		for k in range(self.Z):
			for j in range(self.Y):
				for i in range(self.X):
					count += self.Image[k][j][i]
		VF = float(count)/float(self.Z * self.Y * self.X)
		return VF
	
	def getSize(self):
		""" getSize() is a public method that could be used for every microstructures.
			This function displays the dimension of the microstructure
			@para: None
			@return: None. But print a message of the microstructure size
			Developed on Jul. 15, 2015
		"""
		print " The dimension of the given microstructure is " + str(self.X) + " x " + str(self.Y) + " x " + str(self.Z)
		
	def __str__(self):
		""" default function that could print the microstructure for debugging
			Developed on Jul. 14, 2015
		"""
		print self.Image
		return ""
	
	
	def Coordinates(self, scale):
		""" Coordinates() is a public method that is shared by every microstructure.
			This function could generate the node/element/surface list of the grid that could be fed into FEA
			@para: scale --> the actual length (in micron) of a voxel
			@return: None. But list of nodes/elements/surfaces are generated
			Developed on Jul. 12, 2015
			"""
		Image = self.Image
		Z = self.Z
		Y = self.Y
		X = self.X
		print "== Generating the node list =="
		self.NodeList = []
		NodeIndex = 1
		for k in range(Z+1):
			for j in range(Y+1):
				for i in range(X+1):
					self.NodeList.append([NodeIndex, round(i*scale, 3), round(j*scale, 3), round(k*scale, 3)])
					NodeIndex += 1
		#print self.NodeList
		print "== Generating the element list =="
		self.EleList = []
		Image_Matrix = []
		EleIndex = 1;
		for k in range(Z):
			for j in range(Y):
				for i in range(X):
					self.EleList.append([EleIndex, (k+1)*(X+1)*(Y+1) + (j+1)*(X+1) + i + 1,
						(k+1)*(X+1)*(Y+1)+ j *(X+1) + i + 1,
						k*(X+1)*(Y+1) + j*(X+1) + i + 1,
						k*(X+1)*(Y+1)+(j+1)*(X+1) + i + 1,
						(k+1)*(X+1)*(Y+1)+(j+1)*(X+1) + i + 2,
						(k+1)*(X+1)*(Y+1) + j*(X+1)+ i + 2,
						k*(X+1)*(Y+1) + j*(X+1)+ i + 2,
						k*(X+1)*(Y+1)+(j+1)*(X+1)+ i + 2])
					EleIndex += 1
		print len(self.EleList)
		print "== generating boundary node sets =="
		self.LeftSurface = []
		self.RightSurface = []
		self.TopSurface = []
		self.BotmSurface = []
		self.FrontSurface = []
		self.BackSurface = []
		for k in range(1, Z+2):
			for j in range(1, Y+2):
				self.LeftSurface.append((k-1)*(X+1)*(Y+1)+(j-1)*(X+1)+1)

		for k in range(1, Z+2):
			for j in range(1, Y+2):
				self.RightSurface.append((k-1)*(X+1)*(Y+1)+j*(X+1))

		self.BotmLeftBack = Y*(X+1) + 1
		self.TopLeftBack = Z*(X+1)*(Y+1) + Y*(X+1) + 1
		self.BotmLeftFront = 1
		self.TopLeftFront = Z*(X+1)*(Y+1) + 1

		for j in range(1, Y+2):
			for i in range(1, X+2):
				Vp = Z*(X+1)*(Y+1)+(j-1)*(X+1)+i
				if Vp != self.TopLeftFront:
					self.TopSurface.append(Vp)

		for j in range(1, Y+2):
			for i in range(1, X+2):
				Vp = (j-1)*(X+1) + i
				if Vp != self.BotmLeftFront:
					self.BotmSurface.append(Vp)

		for k in range(1, Z+2):
			for i in range(1, X+2):
				Vp = (k-1)*(X+1)*(Y+1) + i
				if Vp != self.BotmLeftFront:
					self.FrontSurface.append(Vp)

		for k in range(1, Z+2):
			for i in range(1, X+2):
				Vp = (k-1)*(X+1)*(Y+1) + Y*(X+1) + i
				if Vp != self.BotmLeftBack:
					self.BackSurface.append(Vp)
		
		self.scale = scale
		self.Disp = 0.005*scale * self.X
		self.ScaledDisp = self.Disp * scale
		

		print "Node/Surface/Element list generated!"


class individualClusterBase(allClusterBase):
	""" This class is a subclass of allClusterBase. It defines additional method that could be used for individual cluster objects.
		It also inherit the method of its superclass.
	"""

	def add_cluster(self, another):
		""" add_cluster() is a public method that could be used to merge the markers of the same pixel from two individual cluster objects into one.
			@para: another --> has to be a individualClusterBase (or its subclass) object
			@return: None. But self.IntphImage attributes will be updated.
			Developed on Jul. 1st, 2015.
		"""
			
		if self.X != another.X or self.Y != another.Y or self.Z != another.Z: # examine if the dimension matches
			print "Dimension does not match!"
			sys.exit(0)
		else: # if their dimensions match, starts merge structures
			# if a voxel has the impact from two or more clusters, then record the nums of interphases first
			for k in range(self.Z):
				for j in range(self.Y):
					for i in range(self.X):
						if isinstance(self.IntphImage[k][j][i], (int, float)):
							self.IntphImage[k][j][i] = [self.IntphImage[k][j][i], another.IntphImage[k][j][i]]
						else:
							self.IntphImage[k][j][i].append(another.IntphImage[k][j][i])

						# remove duplicate zeros
						# case 1: if this voxel has other label other than 0, remove all zero
						# case 2: if this voxel has multiple zeros, leave only one zero
						zero_count = 0
						pop_list = []
						flag = False
						if 1 in self.IntphImage[k][j][i]: # filler will not be influenced by interphase
							self.IntphImage[k][j][i] = [1]
						for itemIndex in range(len(self.IntphImage[k][j][i])):
							if self.IntphImage[k][j][i][itemIndex] == 0:
								self.IntphImage[k][j][i][itemIndex] = int(self.IntphImage[k][j][i][itemIndex])
								zero_count += 1
								if zero_count > 0:
									pop_list.append(itemIndex)
							else:
								flag = True
						if flag == False:
							pop_list.pop(0)
						pop_list.reverse()
						for index in pop_list:
							self.IntphImage[k][j][i].pop(index)

		print self.IntphImage
	
	def __str__(self):
		""" Temporary toString method for debugging """
		#print self.IntphImage
		return ""
		
	def countIntphLabel(self):
	#countIntphLabel function returns the number of label combinations assigned to the voxels
	
		X = self.X
		Y = self.Y
		Z = self.Z
		max_count = 0
		label_list = []
		for k in range(Z):
			for j in range(Y):
				for i in range(X):
					if isinstance(self.IntphImage[k][j][i], (list,tuple)):
						self.IntphImage[k][j][i].sort()
					label_list.append(self.IntphImage[k][j][i])
		# label_list may have duplicates. 
		NoDup_list = []
		for item in label_list:
			if item not in NoDup_list:
				NoDup_list.append(item) 
		print NoDup_list
		if 0 in NoDup_list:
			#print "There are " + str(len(NoDup_list) - 2) + " Intph overlapping combinations."
			#self.numLabel = len(NoDup_list) - 2
			return	len(NoDup_list) - 2 #exclude 0 and 1, which corresponds to the pure matrix and filler phases
		else:
			#print "There are " + str(len(NoDup_list) - 1) + " Intph overlapping combinations."
			#self.numLabel = len(NoDup_list) - 1
			return len(NoDup_list) - 1
	
	def assignFactors(self, R, BL, BR):
		""" This function is to assign shifting"""
		numLabel = self.countIntphLabel()
		if (len(R) != numLabel) | (len(BL) != numLabel) | (len(BR) != numLabel):
			print "Input array dimension mismatch!"
			sys.exit(0)
		else:
			self.R = R
			self.BL = BL
			self.BR = BR
	
	def assignMatPara(self, Rou_mat, E0_mat, v0_mat, Rou_int, E0_int, v0_int, Rou_Ptc, E_Ptc, v_Ptc, file):
		""" assign materials' parameters"""
		self.Rou_mat = Rou_mat
		self.E0_mat = E0_mat
		self.v0_mat = v0_mat
		self.Rou_mat = Rou_mat
		self.E0_int = E0_int
		self.v0_int = v0_int
		self.Rou_Ptc = Rou_Ptc
		self.E_Ptc = E_Ptc
		self.v_Ptc = v_Ptc
		self.MatFile = file
	
	def assignStepPara(self, f_initial, f_final, NoF, Bias, NumCPU):
		self.f_initial = f_initial
		self.f_final = f_final
		self.NoF = NoF
		self.Bias = Bias
		self.NumCPU = NumCPU
		
	def assignIntphMat(self):
		""" assign Interphase properties based on the processed final label"""
		self.ElePtc = []
		self.EleMat = []
		numLabel = self.countIntphLabel()
		self.EleInt = {}
		for index in range(2, numLabel+2):
			self.EleInt[index] = []
		#print self.EleInt
		Z = self.Z
		Y = self.Y
		X = self.X #make a copy to prevent the key info from accidental changes 
		for k in range(Z):
			for j in range(Y):
				for i in range(X):
					assert isinstance(self.IntphImage[k][j][i], int) # check type
					if self.IntphImage[k][j][i] == 0:
						self.EleMat.append(k*X*Y + j*X + i + 1)
					elif self.IntphImage[k][j][i] == 1:
						self.ElePtc.append(k*X*Y + j*X + i + 1)
					else:
						#print self.EleInt[self.IntphImage[k][j][i]]
						self.EleInt[self.IntphImage[k][j][i]].append(k*X*Y + j*X + i + 1)
		
	
	def composeInput(self, InputName):
		InputFile = open(InputName+".inp", "w")
		self.InputName = InputName
		InputFile.write("*Heading\n")
		InputFile.write("**3D Modeling\n") # description of this model
		InputFile.write('** ----------------------------------------------------------------------- input\n')
		InputFile.write('** PARTS\n')
		InputFile.write('** ----------------------------------------------------------------------- node list\n')
		InputFile.write('*Node\n')
		for i in self.NodeList:
			InputFile.write(str(i)[1:-1]+'\n')
		InputFile.write('**\n')
		InputFile.write('** ----------------------------------------------------------------------- element list\n')
		InputFile.write('*Element, type=C3D8R\n')
		for i in self.EleList:
			InputFile.write(str(i)[1:-1]+'\n')
		InputFile.write('**\n')
		InputFile.write('**------------------------------------------------------------------------ sets\n')
		
		### set boundary node sets           
		InputFile.write('*Nset, nset=LeftSurface\n')
		for i in self.LeftSurface:
			InputFile.write(str(i)+'\n')
			
		InputFile.write('*Nset, nset=RightSurface\n')
		for i in self.RightSurface:
			InputFile.write(str(i)+'\n')
			
		InputFile.write('*Nset, nset=TopSurface\n')
		for i in self.TopSurface:
			InputFile.write(str(i)+'\n')
			
		InputFile.write('*Nset, nset=BottomSurface\n')
		for i in self.BotmSurface:
			InputFile.write(str(i)+'\n')
			
		InputFile.write('*Nset, nset=FrontSurface\n')
		for i in self.FrontSurface:
			InputFile.write(str(i)+'\n')
			
		InputFile.write('*Nset, nset=BackSurface\n')
		for i in self.BackSurface:
			InputFile.write(str(i)+'\n')

		InputFile.write('*Nset, nset=BotmLeftBack\n')  
		InputFile.write(str(self.BotmLeftBack)+'\n')

		InputFile.write('*Nset, nset=TopLeftBack\n')  
		InputFile.write(str(self.TopLeftBack)+'\n')

		InputFile.write('*Nset, nset=BotmLeftFront\n')  
		InputFile.write(str(self.BotmLeftFront)+'\n')

		InputFile.write('*Nset, nset=TopLeftFront\n')  
		InputFile.write(str(self.TopLeftFront)+'\n')
		
		InputFile.write('*Elset, elset=Set-PTC\n') 
		for i in self.ElePtc:
			InputFile.write(str(i)+'\n')
		
		numLabel = self.countIntphLabel()
		for index in range(2, numLabel+2):
			InputFile.write('*Elset, elset=Set-INT-' + str(index-1) + '\n') 
			for i in self.EleInt[index]:
				InputFile.write(str(i)+'\n')

		#InputFile.write('*Elset, elset=Set-INT-2\n') 
		#for i in self.EleInt2:
		#	InputFile.write(str(i)+'\n')

		#InputFile.write('*Elset, elset=Set-INT-3\n') 
		#for i in self.EleInt3:
		#	InputFile.write(str(i)+'\n')
		
		#InputFile.write('*Elset, elset=Set-MAT\n')
		#for i in self.EleMat:
		#	InputFile.write(str(i)+'\n')
		InputFile.write('*Elset, elset=Set-MAT\n')
		for i in self.EleMat:
			InputFile.write(str(i)+'\n')
		
		InputFile.write('** -------------------------------------------------------------------- section\n') 
		InputFile.write('** Section: Sec-MAT\n')
		InputFile.write('*Solid Section, elset=Set-MAT, material=MAT-MATRIX\n')
		InputFile.write('1.,\n')
		for index in range(1, numLabel+1):
			InputFile.write('** Section: Sec-INT-'+ str(index) +'\n')
			InputFile.write('*Solid Section, elset=Set-INT-' +str(index) +', material=MAT-INTERPHASE-1\n')  ## first layer of interphase
			InputFile.write('1.,\n')
#InputFile.write('** Section: Sec-INT-2\n')
#InputFile.write('*Solid Section, elset=Set-INT-2, material=MAT-INTERPHASE-2\n')  ## second layer of interphase
#InputFile.write('1.,\n')
#InputFile.write('** Section: Sec-INT-3\n')
#InputFile.write('*Solid Section, elset=Set-INT-3, material=MAT-INTERPHASE-3\n')  ## third layer of interphase
#InputFile.write('1.,\n')
		InputFile.write('** Section: Sec-PTC\n')
		InputFile.write('*Solid Section, elset=Set-PTC, material=MAT-PARTICLE\n')
		InputFile.write('1.,\n')
		
		InputFile.write('** --------------------------------------------------------------------- equation\n')
		### set constraints
		InputFile.write('** Constraint: Eqn-1\n')
		InputFile.write('*Equation\n')
		InputFile.write('2\n')
		InputFile.write('TopSurface, 3, 1.\n')
		InputFile.write('TopLeftFront, 3, -1.\n')
		InputFile.write('** Constraint: Eqn-2\n')
		InputFile.write('*Equation\n')
		InputFile.write('2\n')
		InputFile.write('BottomSurface, 3, 1.\n')
		InputFile.write('BotmLeftFront, 3, -1.\n')
		InputFile.write('** Constraint: Eqn-3\n')
		InputFile.write('*Equation\n')
		InputFile.write('2\n')
		InputFile.write('FrontSurface, 2, 1.\n')
		InputFile.write('BotmLeftFront, 2, -1.\n')
		InputFile.write('** Constraint: Eqn-4\n')
		InputFile.write('*Equation\n')
		InputFile.write('2\n')
		InputFile.write('BackSurface, 2, 1.\n')
		InputFile.write('BotmLeftBack, 2, -1.\n')
		
		### set materials
		InputFile.write('**---------------------------------------------------------------------- material\n') 
		InputFile.write('**\n')
		InputFile.write('** MATERIALS\n')
		InputFile.write('**\n')
		
		InputFile.write('*Material, name=MAT-PARTICLE\n')
		InputFile.write('*Density\n')
		InputFile.write(' '+str(self.Rou_Ptc)+',\n') 
		InputFile.write('*Elastic\n')
		InputFile.write(' '+str(self.E_Ptc)+', '+str(self.v_Ptc)+'\n')
		
		MatrlFile = open(self.MatFile)
		lines = MatrlFile.readlines()
		freq = []
		Ep = []
		Epp = []
		TanDelta = []
		for ln in lines[0:len(lines)]:
			line=ln.strip().split()
			freq.append(float(line[0]))
			Ep.append(float(line[1]))
			Epp.append(float(line[2]))
			TanDelta.append(float(line[2])/float(line[1]))
		freq.reverse()
		Ep.reverse()
		Epp.reverse()
		TanDelta.reverse()
		m = max(TanDelta)
		location = TanDelta.index(m)
		E0 = max(Ep)  ## Now transform Young's Modulus into Shear Modulus assuming Bulk Modulus is constant
		K0 = E0/(3*(1-2*self.v0_mat))
		E_inft = min(Ep) ## long time Young's Modulus
		v_inft = (3*K0-E_inft)/(6*K0)
		Gp = []  ## G'
		Gpp = [] ## G"
		for i in range(len(Ep)):
			EComp = complex(Ep[i],Epp[i])
			GComp = 3*K0*EComp/(9*K0-EComp)
			Gp.append(GComp.real)
			Gpp.append(GComp.imag)
		G_inft = min(Gp)
		
		InputFile.write('*Material, name=MAT-MATRIX\n')
		InputFile.write('*Density\n')
		InputFile.write(' '+str(self.Rou_mat)+',\n')
		InputFile.write('*Elastic, moduli=LONG TERM\n')
		InputFile.write(' '+str(E_inft/10**12)+', '+str(v_inft)+'\n')
		InputFile.write('*Viscoelastic, frequency=TABULAR\n')
		for i in range(len(Gp)):
			InputFile.write(str(Gpp[i]/G_inft)+',\t' +str(1-Gp[i]/G_inft)+',\t' +'0,\t'+'0,\t' + str(freq[i]/(2*math.pi)) + '\n')
		
		for index in range(1, numLabel+1):
		# Interphase properties
			InputFile.write('*Material, name=MAT-INTERPHASE-'+str(index)+'\n') ## the #index layer of interphase 
			InputFile.write('*Density\n')
			InputFile.write(' '+str(self.Rou_mat)+',\n')
			InputFile.write('*Elastic, moduli=LONG TERM\n')
			InputFile.write(' '+str(E_inft/10**12)+', '+str(v_inft)+'\n')
			InputFile.write('*Viscoelastic, frequency=TABULAR\n')
			f = [math.log10(i)-self.R[index-1] for i in freq]
			
			fb=[]
			for i in range(location):
				f_b=(f[i]-f[location])*self.BL[index-1]+f[location]
				fb.append(f_b)
			for i in range(location,len(Gp)):
				f_b=(f[i]-f[location])*self.BR[index-1]+f[location]
				fb.append(f_b)
			newfreq = [10**i for i in fb]
			
			for i in range(len(Gp)):
				InputFile.write(str(Gpp[i]/G_inft)+',\t' +str(1-Gp[i]/G_inft)+',\t' +'0,\t'+'0,\t' + str(newfreq[i]/(2*math.pi)) + '\n')
		
		MatrlFile.close() # close the master curve file
		
		## set loading step
		InputFile.write('** ----------------------------------------------------------------\n')
		InputFile.write('*STEP,NAME=STEP-1,PERTURBATION\n')
		InputFile.write('3D,TENSIONXX\n')
		InputFile.write('*STEADY STATE DYNAMICS, DIRECT\n')
		InputFile.write(str(self.f_initial) + ', ' + str(self.f_final) + ', ' + str(self.NoF) + ', ' + str(self.Bias) + ',\n')
		InputFile.write('*BOUNDARY\n')
		InputFile.write('BotmLeftFront,1,6\n')
		InputFile.write('LeftSurface, 1, 1\n')
		#InputFile.write('LeftSurface,1,6\n')
		#InputFile.write('RightSurface, 2, 6\n')
		InputFile.write('RightSurface, 1, 1, ' + str(self.Disp) + '\n')
		InputFile.write('***************************************************************************************************\n')
		InputFile.write('**OUTPUT REQUEST\n')
		InputFile.write('**\n')
		InputFile.write('*RESTART,WRITE,frequency=0\n')
		InputFile.write('**\n')
		InputFile.write('**FIELD OUTPUT: F-OUTPUT-1\n')
		InputFile.write('**\n')
		InputFile.write('*Output, field, variable=PRESELECT\n')
		InputFile.write('*element output\n')
		InputFile.write('evol\n')
		InputFile.write('**\n')
		InputFile.write('**HISTORY OUTPUT: H-OUTPUT-1\n')
		InputFile.write('**\n')
		InputFile.write('*END STEP\n') 
		InputFile.close()
		
		print 'Input file has been successfully generated'
		
		
	def runJob(self):
		""" By using this function, the .inp file will be fed into ABAQUS to run 
		@para: None
		@return: None
		Developed on Jul. 18, 2015"""
		assert hasattr(self, 'InputName')
		assert hasattr(self, 'NumCPU')
		assert os.path.isfile(self.InputName + ".inp")
		print "Now submitting the input file into ABAQUS for analysis..."
		os.system('abq6112 job='+self.InputName+' cpus='+str(self.NumCPU)+ ' int')
		
		print "Job completed!"
		
	def extractData(self):
		assert hasattr(self, 'InputName')
		readODB = open('readODB'+str(self.InputName)+'.py','w')
		readODB.write('import math\n')
		readODB.write('from odbAccess import *\n')
		readODB.write('from abaqusConstants import *\n')
		readODB.write('InputName=\''+str(self.InputName)+'\'\n')
		readODB.write('print \'ODB = '+str(self.InputName) +'\'\n')
		readODB.write('outputfilename=InputName+\'-Modulii.txt\'\n')
		readODB.write('out = open(outputfilename,\'w\')\n')
		readODB.write('odb = openOdb(InputName+\'.odb\')\n')
		readODB.write('numberFrame = len(odb.steps[\'STEP-1\'].frames)\n')
		readODB.write('nset=odb.rootAssembly.instances[\'PART-1-1\'].nodeSets[\'RIGHTSURFACE\']\n')
		readODB.write('for ns in range (1,numberFrame):\n')
		readODB.write('\tframe = odb.steps[\'STEP-1\'].frames[ns]\n')
		readODB.write('\tforce=frame.fieldOutputs[\'RF\']\n')                                                	
		readODB.write('\tnsetForce = force.getSubset(region=nset)\n')                                      	
		readODB.write('\tsumforceReal=0\n')                                                                	
		readODB.write('\tsumforceImag=0\n')                                                                	
		readODB.write('\tfor v in nsetForce.values:\n')
		readODB.write('\t\tsumforceReal+=v.data[0]/('+str(self.ScaledDisp)+'*'+str(self.Y)+'*'+str(self.Z)+'/'+str(self.X)+')\n')
		readODB.write('\t\tsumforceImag+=v.conjugateData[0]/('+str(self.ScaledDisp)+'*'+str(self.Y)+'*'+str(self.Z)+'/'+str(self.X)+')\n')                             	
		readODB.write('\tout.write(str(frame.frameValue*2*math.pi)+\'\\t\'+str(sumforceReal)+\'\\t\'+str(sumforceImag)+\'\\n\')'+'\n')
		readODB.write('out.close()\n')
		readODB.write('odb.close()\n')
		readODB.close()
		print 'Now retrieving modulus from ODB file'
		os.system('abq6112 python readODB'+str(self.InputName)+'.py')
		print 'Modulus file generated !!'
		
		
	#==========================================================
	#=           Python Interfaces (Abstract method)          =
	#==========================================================
	def mergeIntph(self):
		""" Interface that to be implemented to merge the markers in the pixels' Intph fiels into one.
			Has to be implemented in the subclass to avoid error.
			Developed on Jul 15, 2015.
		"""
		raise NotImplementedError()
		
	def assignIntph(self, num_intph):
		raise NotImplementedError()
		
		

class clusterMethod1(individualClusterBase):
	""" this class is the subclass of individualClusterBase that provides implementations of the two interfaces (mergeIntph, assignIntph)
	"""
	
	def assignIntph(self, num_intph):
		""" assignIntph() is the implementation of the interface assignIntph() in the base function.
			It generates num_intph intph layers around the cluster by extending the cluster boarder.
			@para: num_intph --> the number of interphase layers that we want to add
			@return: None, but will update the self.IntphImage field.
			Developed on Jul. 16, 2015
		"""
		
		X = self.X
		Y = self.Y
		Z = self.Z
		Image = numpy.array(self.Image, dtype=object)
		# find the current outer most layer
		max_counter = 0;
		for z in range(Z):
			for y in range(Y):
				for x in range(X):
					if Image[z][y][x] > max_counter:
						max_counter = Image[z][y][x]
		
		# add interphase layers
		current_layer = max_counter
		new_layer = max_counter + 1;
		while current_layer < num_intph + max_counter:
			for z in range(1, Z+1):
				for y in range(1, Y+1):
					for x in range(1, X+1):
						if Image[z-1][y-1][x-1] == current_layer: 
						#find a outer most cluster pixel, update the neighboring matrix(0) pixel with the new layer label
							for k in range(-1, 2):
								for j in range(-1, 2):
									for i in range(-1, 2):
										if x+i>=1 and x+i<=X and y+j>=1 and y+j<=Y and z+k>=1 and z+k<=Z:
											if Image[z-1+k][y-1+j][x-1+i] == 0:
												Image[z-1+k][y-1+j][x-1+i] = int(new_layer)
			current_layer += 1
			new_layer += 1
		print Image
		self.IntphImage = Image
		
		
	def mergeIntph(self):
		""" mergeIntph() is an implementation of the interface in the superclass 
			This function selects the minimun marker in the list which indicates the nearest-neighbor-surface-distance assumption
			Developed on Jul. 16, 2015
			"""
		X = self.X
		Y = self.Y
		Z = self.Z
		for z in range(Z):
			for y in range(Y):
				for x in range(X):
					self.IntphImage[z][y][x] = min([item for item in self.IntphImage[z][y][x]])
		#print self.IntphImage


