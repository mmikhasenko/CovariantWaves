# original author Fabian Krinner (2016)


import sympy
from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorType, tensorsymmetry, tensorhead
from sympy.functions.special.tensor_functions import eval_levicivita
from sympy.physics.quantum.cg import CG
from sympy import symbols, simplify, I, Rational, factorial
import itertools

Lorentz = TensorIndexType('Lorentz', dummy_fmt='L', dim = 4, eps_dim = 4)
Lorentz.data = [1,-1,-1,-1]
epsilonData = [[[[0 for _ in range(4)] for __ in range(4)] for ___ in range(4)] for ____ in range(4)]
for i in range(4):
	for j in range(4):
		for k in range(4):
			for l in range(4):
				epsilonData[i][j][k][l] = eval_levicivita(i,j,k,l)

Lorentz.epsilon.data = epsilonData

symV = tensorsymmetry([1])
vType = TensorType([Lorentz], symV)
symT = tensorsymmetry([1,1])
tType = TensorType([Lorentz, Lorentz], symT)

beta, gamma = symbols("beta gamma")

E0Rest = vType("E0")
E0Rest.data = [0,0,0,1]
#E0Rest.data = [beta*gamma,0,0,gamma]
EpRest = vType("Ep")
EpRest.data = [0,-1/sympy.sqrt(2),-I/sympy.sqrt(2),0]
EmRest = vType("Em")
EmRest.data = [0,1/sympy.sqrt(2),-I/sympy.sqrt(2),0]

extractors = []
for i in range(4):
	adder = vType("e_"+str(i))
	data = [0,0,0,0]
	if i == 0:
		data[i] = 1
	else:
		data[i] =-1
	adder.data = data
	extractors.append(adder)

globalIndexCounts = {}
def getNewLorentzIndex(baseCharacter = 'm'):
	global globalIndexCounts
	if not baseCharacter in globalIndexCounts:
		globalIndexCounts[baseCharacter] = 0


	i = tensor_indices(baseCharacter+'_'+str(globalIndexCounts[baseCharacter]), Lorentz)
	globalIndexCounts[baseCharacter] += 1
	return i

def getX(kp, gp, J):
	indices = [getNewLorentzIndex() for _ in range(J)]
	dim = J
	if dim == 0:
		return 1
	if dim == 1:
		return kp(indices[0])
	else:
		first = True
		for i in range(dim):
			indexSubList = []
			for ii in range(dim):
				if not i == ii:
					indexSubList.append(indices[ii])
			lessTens = getX(kp, gp, J-1)
			addTens = kp(indices[i])*lessTens(*indexSubList)
			if first:
				retTens1 = addTens
				first = False
			else:
				retTens1 += addTens
		first = True
		for i in range(0,dim-1):
			for j in range(i+1, dim):
				indexSubList = []
				for ii in range(dim):
					if not i == ii and not j == ii:
						indexSubList.append(indices[ii])
				lessTens = getX(kp, gp, J-2)
				if dim == 2:
					addTens = gp(indices[i], indices[j])
				else:
					addTens = gp(indices[i], indices[j])*lessTens(*indexSubList)
				if first:
					retTens2 = addTens
					first = False
				else:
					retTens2 += addTens
	iTens = getNewLorentzIndex()
	kp2 = kp(iTens)*kp(-iTens)
	retTens = ((2*dim-1) * retTens1 - 2*kp2 * retTens2)/dim/dim
	return simplify(retTens)

def getPolarizatrionVector(M, P):
	symV = tensorsymmetry([1])
	vType = TensorType([Lorentz], symV)
	if M == 1:
		return EpRest
	if M == -1:
		return EmRest
	if M == 0:
		return E0Rest
	else:
		raise ValueError("M = "+str(M)+ " not possible for vector")
	return polVec

def getPolarizatrionTensorRest(J, M, P = 22): #P is not needed
	if M > J or M < -J:
		raise ValueError("spin "+str(J)+" state can't have M = "+str(M))
	if J == 0:
		return 1
	if J == 1:
		return getPolarizatrionVector(M,P)
	elif J > 1:
		j = J - 1
		jIndices = [getNewLorentzIndex() for _ in range(j)]
		onedex = getNewLorentzIndex()
		first = True
		for m1 in [-1,0,1]:
			for m2 in range(-j,j+1):
				if m1 + m2 == M:
					toAdd = CG(1,m1,j,m2,J,M).doit() * getPolarizatrionTensorRest(j,m2,P)(*jIndices) * getPolarizatrionVector(m1,P)(onedex)
					if first:
						first  = False
						retVal = toAdd
					else:
						retVal+= toAdd
		return retVal

def ars(s,r):
	if s < 2*r:
		raise ValueError(" s < 2*r in ars")
	val = Rational(-1,2)**r*factorial(s)/factorial(r)/factorial(s-2*r)
	for i in range(1,r+1):
		val /= 2*s - 2*i + 1
#	print "ars gives:",val
	return val


def get_gabgab(gp, aIndices, bIndices):
	if not len(aIndices) == len(bIndices):
		raise Exception("number of a and b indices do not match")
	if len(aIndices) == 0:
		return 1
	retVal = gp(aIndices[0], bIndices[0])
	for i in range(1, len(aIndices)):
		retVal *= gp(aIndices[i], bIndices[i])
	return retVal

def get_gaagbb(gp, aIndices, bIndices):
	if not len(aIndices) == len(bIndices):
		raise Exception("number of a and b indices do not match")
	if len(aIndices) == 0:
		return 1
	if not len(aIndices)%2 == 0:
		raise Exception("get_gaagbb needs an even number of indices")
	retVal = gp(aIndices[0], aIndices[1])*gp(bIndices[0], bIndices[1])
	for i in range(1, len(aIndices)//2):
		retVal *= gp(aIndices[2*i], aIndices[2*i+1])*gp(bIndices[2*i], bIndices[2*i+1])
	return retVal

def getUnymmetrizedProjector(gp, indices):
	if not len(indices)%2 == 0:
		raise Exception("odd number of indices...")
	S = len(indices)//2
	aIndices = indices[S:]
	bIndices = indices[:S]
	if not len(aIndices) == len(bIndices):
		raise Exception("number of a and b indices do not match")
	retVal = get_gabgab(gp, aIndices, bIndices)

	for r in range(1,S/2+1):
		retVal += ars(S, r)*get_gaagbb(gp, aIndices[:2*r], bIndices[:2*r])*get_gabgab(gp,aIndices[2*r:], bIndices[2*r:])
	return retVal

def getProjector(S, gp):
	aIndices = [getNewLorentzIndex('a') for a in range(S)]
	bIndices = [getNewLorentzIndex('b') for b in range(S)]
	unSymm = getUnymmetrizedProjector(gp, aIndices+bIndices)
	first = True
	aPerms = list(itertools.permutations(aIndices, S))
	bPerms = list(itertools.permutations(bIndices, S))
	for a,aPerm in enumerate(aPerms):
		print('a',a,len(aPerms))
		aaa = [aPerm[i] for i in range(S)]
		totalIndices = aaa + bIndices
		if first:
			intermed = unSymm(*totalIndices)
			first = False
		else:
			intermed += unSymm(*totalIndices)
	first = True
	for b,bPerm in enumerate(bPerms):
		print('b',b,len(bPerms))
		bbb = [bPerm[i] for i in range(S)]
		totalIndices = aIndices + bbb
		if first:
			retVal = intermed(*totalIndices)
			first = False
		else:
			retVal += intermed(*totalIndices)
	retVal /= factorial(S)**2
	return retVal

def contract(tens1, tens2, j, P):
	try:
		i1 = tens1.rank
	except AttributeError:
		i1 = 0
	try:
		i2 = tens2.rank
	except AttributeError:
		i2 = 0
	if j > i1+i2 or j < abs(i1-i2):
		return 0
	if (i1+i2-j)%2==0:
		k = (i1+i2-j)//2
		indices1 = [getNewLorentzIndex() for _ in range(i1 - k)]
		indices2 = [getNewLorentzIndex() for _ in range(i2 - k)]
		for _ in range(k):
			iContr = getNewLorentzIndex()
			indices1.append( iContr)
			indices2.append(-iContr)
		if i1==0:
			if i2==0:
				return tens1 * tens2
			else:
				return tens1*tens2(*indices2)
		else:
			return tens1(*indices1)*tens2(*indices2)
	else:
		k = (i1+i2-j-1)//2
		indices1 = []
		indices2 = []
		for _ in range(k):
			iContr = getNewLorentzIndex()
			indices1.append( iContr)
			indices2.append(-iContr)
		epsilonIndices = [getNewLorentzIndex() for _ in range(4)]
		indices1.append(-epsilonIndices[0])
		indices2.append(-epsilonIndices[1])
		while len(indices1) < i1:
			indices1.append(getNewLorentzIndex())
		while len(indices2) < i2:
			indices2.append(getNewLorentzIndex())
		return tens1(*indices1)*tens2(*indices2)*P(-epsilonIndices[2])*Lorentz.epsilon(*epsilonIndices)
	
def fullContract(tens1, tens2):
	try:
		i1 = tens1.rank
	except AttributeError:
		i1 = 0
	try:
		i2 = tens2.rank
	except AttributeError:
		i2 = 0
	if not i1 == i2:
		raise Exception("ranks do not match")
	if i1 == 0:
		return tens1*tens2
	ind1 = []
	ind2 = []
	for i in range(i1):
		index = getNewLorentzIndex()
		ind1.append( index)
		ind2.append(-index)
	return tens1(*ind1)*tens2(*ind2)


def extract_values(tensor):
	if isinstance(tensor, int):
		return tensor
	dim = tensor.rank
	retter = []
	indices = [getNewLorentzIndex() for _ in range(dim)]
	for i in range(4):
		if dim > 1:
			retter.append(extract_values(extractors[i](-indices[0])*tensor(*indices)))
		else:
			retter.append((extractors[i](-indices[0])*tensor(*indices)).data)
	return retter

def makeCstring(J ,M, P):
	tens = getPolarizatrionTensorRest(J,M,P)
	string = str(extract_values(tens)).replace('[','{').replace(']','}')
	for i in range(J+1):
		expStringGamma = '*'.join(['gamma']*i)
		expStringBeta  = '*'.join(['beta']*i)
		string = string.replace('gamma**'+str(i), expStringGamma).replace('beta**'+str(i), expStringBeta)

	base = "std::complex<double> phi"+str(J)+str(M).replace('-','m')+" "
	for i in range(J):
		base+='[4]'
	base += ' = '+string+';'
	return base

def main():
	i, j, k, l = tensor_indices('i,j,k,l', Lorentz)
	kp = vType('kp')
	kp0,kp1,kp2,kp3 = symbols("kp0 kp1 kp2 kp3")
	kp.data = [kp0,kp1,kp2,kp3]

	P = vType('P')
#	p0,p1,p2,p3 = symbols("p0 p1 p2 p3")
#	P.data = [p0,p1,p2,p3]
#	P2 = symbols("M3pi2")

#	gp = tType("gp")
#	gp.data = [[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]]
#	gp = gp(i,j) - P(i)*P(j)/(P**2)

	kp = vType("kp")
	gp = tType("gp")
	
	X3 = getX(kp,gp,3)
	simplify(X3)
	print(X3)

	gp = tType('gp')
#	X1 = getX(kp,gp,3)
#	X2 = getX(kp,gp,2)
#	J = 2
#	TT = contract(X1, X2, J, P)
#	pol = getPolarizatrionTensor(J,1,P)

#	with open("theCprogram.cc", 'w') as out:
#		out.write("#include<complex>\nint main(){\n\tdouble beta;\n\tdouble gamma;\n\tstd::complex<double> I(0.,1.);\n")
#		for J in range(7):
#			print "start J =",J
#			for M in range(-J,J+1):
#				print "M =",M
#				out.write('\t'+makeCstring(J,M,P)+'\n')
#		out.write("\treturn 0;\n}")

#	print getProjector(1,gp)

if __name__ == "__main__":
	main()



