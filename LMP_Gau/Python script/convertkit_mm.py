import sys
import os
import numpy as np
import lammps
import openbabel as ob



class convertkit:
		def __init__(self, inputfile=None):
		# wishlist:
		#	mol.readG09Output(filename='c3h4o5_ts_bla.fchk')
		#	mal testen fur self.X
		#		class C(object):
		#		def __init__(self):
		#			self._x = None
		#
		#		def _x_get(self):
		#			return self._x
		#
		#		def _x_set(self, value):
		#			self._x = value
		#
		#		def _x_del(self):
		#			del self._x
		#
		#		x = property(_x_get, _x_set, _x_del, "I'm the 'x' property.")
		#   orientate() -> ursprung in com, xyz auf Ip
		#   self.w etc initialisieren als array([0,0,0])
		#	find suitable convergence criterium for IRC
		#	simplify hessian calculation
		#
			# . files
			#self.ffield = ffield
			#self.title = title if title else os.path.basename(inputfile) if inputfile else 'notitle'
			self.logflag = False
			self.evoutflag = False # True means file has (already) been written to
			self.ircoutflag = False
			self.bondfilename = 'tmp.bonds'

			# . gaussian stuff
			self.qmcount = 0
			self.spawncount = 0

			# . openbabel converter (for smiles)
			ob.obErrorLog.SetOutputLevel(0)
			self.conv = ob.OBConversion()

			# . fixed molecule info
			self.type_to_number = {'C':6,'H':1,'O':8,'N':7}
			self.mass = {1: 1.007825032, 6: 12.0, 7: 14.0, 8: 15.999}
			self.types = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
			self.valence = {1: 1, 2: 0, 6: 4, 7: 5, 8: 2, 9: 1, 10: 0, 16: 6, 17: 1, 18: 0} # z -> valence
			self.coordination = {1: 1, 2: 0, 6: 4, 7: 3, 8: 2, 9: 1, 10: 0, 16: 6, 17: 1, 18: 0}
			self.a = 50 # irgendwas > 10
			self.atomtypes = 0
			self.identify = {}
			#self.active = active if active!=None else []
			self.atomicnum = []

			# . molecule state
			self.X = np.array([])
			self.X_was_reset = True
			self.H = False
			self.G = False
			self.e = False
			self.w = False
			self.V = False
			#self.k_follow = k_follow
			self.rmsf = False
			self.maxf = False
			self.converged = False
			self.graph     = False
			self.molecules = False
			self.mom       = False
			self.smiles    = False
			self.Q = 1
			# . the 3 exes
			#   TS, Min1 and Min2 coordinates
			self.Xinit = np.array([])
			self.Xts   = np.array([])
			self.Xmin1 = np.array([])
			self.Xmin2 = np.array([])

			# . algorithm options
			self.Slimit = 0.01
			self.omitev = True
			self.Wlimit = 0.1 # 0.5
			self.follow_min = False
			#self.hardopt = False
			
			# . hessian accuracy. minimum 1e-6
			self.hessian_dx = 1e-5
			
			# . Trust Radius Method
			self.Tradius = 0.03 # 0.01588
			self.Tlimit  = 0.1
			self.alpha   = 1.5
			self.q_good  = 0.85
			self.q_min   = 0.70
			self.Ntruststeps = 10

			# . default convergence limits
			#   from gaussian09 'tight' convergence:
			#    F    = 0.0003  Ha/Bohr * 1185.82  = 0.356 kcal/mole-A
			#    Fmax = 0.00045 Ha/Bohr * 1185.82  = 0.534 kcal/mole-A
			#    S    = 0.0012  Bohr    * 0.529177 = 0.000635 A
			#    Smax = 0.0018  Bohr    * 0.529177 = 0.000953 A
			#   but looser limit for steps, since algorithm dependent
			self.Fconv    = 0.0003  * 1185.82  # = 0.356 kcal/mole-A
			self.Fmaxconv = 0.00045 * 1185.82  # = 0.534 kcal/mole-A
			self.Sconv    = 0.12    * 0.529177 # = 6.35e-1 A
			self.Smaxconv = 0.18    * 0.529177 # = 9.53e-1 A

		# . load XYZ
			if inputfile != None:
				self.load(inputfile)
			
			
			
		def load(self, Infile, Write=True):
				    reader = open(Infile)
				    line = reader.readline()
				    if line.startswith('%'):
					    reader.close()
					    self.loadCOM(Infile,Write=Write)
				    else:
					    reader.close()
					    self.loadXYZ(Infile,Write=Write)
				    return



		def loadXYZ(self, XYZfile, Write=False):
			  print '===', 'load', XYZfile, '==='
			  reader = open(XYZfile)

			  # . read number of atoms
			  self.N = int(reader.readline().strip().split()[0])

			  # . read, print and process comment
			  comment = reader.readline()
			  print 'Comment:', comment.strip()
			  for p in comment.split():
				  if p.startswith('follow='):
					  self.k_follow = int(p.split('=')[1]) - 1
					  print 'set follow EV to', self.k_follow
				  if p.startswith('active='):
					  self.active = [int(i) for i in p.split('=')[1].strip('{}[]()').split(',')]
					  print 'set active atoms to', self.active

			  # . pack geometry
			  geo = []
			  self.atomicnum = []
			  ids = {}
			  for i in xrange(self.N):
				  words = reader.readline().strip().split()
				  z = words[0]
				  # . handle atom type, e.g. '6' or 'C' is possible
				  if z.isdigit():
					  z = int(z)
					  if z not in [1,6,7,8]: print 'WARNING: Atom type in input not C H O or N!'
				  else:
					  if z not in ['C','H','O','N']:
						  print 'WARNING: Atom type in input not C H O or N!'
						  z = 1
					  else:
						  z = self.type_to_number[z]
				  self.atomicnum.append(z)
				  ids[z] = True
				  geo += [float(words[1]), float(words[2]), float(words[3])]

			  # . save mol info
			  self.converged = False
			  self.X = np.array(geo)
			  self.Xinit = self.X.copy()
			  self.atomtypes = len(ids)
			  self.identify = {sorted(ids)[i]: i+1 for i in xrange(self.atomtypes)}

			  # . initialize energy, gradient
			  #self.calcForces()

			  # . print energy info
			  print 'Energy:      %+14f (Kcal/mole)' %(self.e)
			  print 'Force RMS:   %+14f (Kcal/mole-Angstrom)' %(self.rmsf)
			  print 'Max Force:   %+14f (Kcal/mole-Angstrom)' %(self.maxf)

			  # . flag as converged
			  self.converged = ( self.rmsf < self.Fconv and self.maxf < self.Fmaxconv )

			  # . log initial structure
			  #if Write: self.logGeometry(Comment='initial: %.2f kcal/mol %.3f kcal/mol/A' %(self.e,self.rmsf))

			  reader.close()
			  return

		def logGeometry(self, Comment=None, Geo=None):
	# . write an xyz file entry
	#
			if Geo == None: Geo = self.X
			if Comment == None: Comment = '%s: %.2f kcal/mol %.3f kcal/mol/A' %(os.path.basename(self.ffield),self.e,self.rmsf)
			if not self.logflag:
				self.logflag = True
				self.logxyz = open(self.title+'.nr_out.xyz', 'w', 1)
			self.logxyz.write('%d\n' %(self.N,))
			self.logxyz.write('%s\n' %(Comment,))
			for i in xrange(self.N):
				self.logxyz.write('%s %f %f %f\n' %(self.types[self.atomicnum[i]], Geo[3*i+0], Geo[3*i+1], Geo[3*i+2]))
			return

		def logX(self, Comment=None, X=None):
		# . shortcut
		#
			self.logGeometry(Comment,X)
			return

		def writeLammpsDataFile(self, X=None, atomicnum=None, filename='tmp.data'): #where is tmp.data created...? in molkit tmp is implemented in def canonical (tmp = obmol2.GetBond(atom, btom))??
			  if X == None: X = self.X
			  if atomicnum == None: atomicnum = self.atomicnum
			  datafile = open(filename, 'w')
			  datafile.write('\n')
			  datafile.write(' %d atoms\n' %(len(X)/3))
			  datafile.write('\n')
			  datafile.write(' %d atom types\n' %(self.atomtypes))
			  datafile.write('\n')
			  datafile.write(' %.02f %.02f xlo xhi\n' %(-self.a, self.a))
			  datafile.write(' %.02f %.02f ylo yhi\n' %(-self.a, self.a))
			  datafile.write(' %.02f %.02f zlo zhi\n' %(-self.a, self.a))
			  datafile.write('\n')
			  datafile.write(' Masses\n')
			  datafile.write('\n')
			  for key in sorted(self.identify):
				  datafile.write(' %d %f\n' %(self.identify[key], self.mass[key]))
			  datafile.write('\n')
			  datafile.write(' Atoms\n')
			  datafile.write('\n')
			  for i in xrange(len(X)/3):
				  datafile.write(' %d 0 %d %f %f %f\n' %(i+1, self.identify[atomicnum[i]], X[3*i+0], X[3*i+1], X[3*i+2] ) )
			  datafile.write('\n')
			  #datafile.close()
			  print(datafile)
			  return
			
			
			
