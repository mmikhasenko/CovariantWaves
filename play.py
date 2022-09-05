from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorType, tensorsymmetry, tensorhead
from sympy import symbols
Lorentz = TensorIndexType('Lorentz', dummy_fmt='L', dim = 4, eps_dim = 4)
Lorentz.data = [1,-1,-1,-1]

symV = tensorsymmetry([1])
vType = TensorType([Lorentz], symV)
symT = tensorsymmetry([1,1])
tType = TensorType([Lorentz, Lorentz], symT)

globalIndexCounts = {}
def getNewLorentzIndex(baseCharacter = 'm'):
        global globalIndexCounts
        if not baseCharacter in globalIndexCounts:
                globalIndexCounts[baseCharacter] = 0


        i = tensor_indices(baseCharacter+'_'+str(globalIndexCounts[baseCharacter]), Lorentz)
        globalIndexCounts[baseCharacter] += 1
        return i

def isDummy(index):
	name = index.name
	if not name.startswith('dummy_index_'):
		return index.is_up, -1
	return index.is_up, int(name[12:])

def refineIndices(indices, infos):
	nInd = len(indices)
	if not nInd == len(infos):
		raise Exception("size mismatch in indices and indofs")
	already = []
	double = []
	for i in range(nInd):
		ind = infos[i][1]
		if ind in double:
			raise Exception("index found trice (three times)")
		if ind in already:
			double.append(ind)
		already.append(ind)
	if len(double) == 0:
		return indices, infos
	retInd = []
	retInf = []
	for i in range(nInd):
		ind = infos[i][1]
		if not ind in double:
			retInd.append(indices[i])
			retInf.append(infos[i])
	return retInd, retInf
		

def applyIndetityToSingeExpression(identity, expression):
	retVal = expression[0]
	countIndices = 0
	same = identity[0] == identity[1]
	firstList  = []
	secondList = []
	dummyIndexInfos1 = []
	dummyIndexInfos2 = []
	for i in range(len(expression[1])):
		nInd = expression[1][i].rank
		if nInd > 1:
			raise Exception("contraction of rank 2 objects not implemented yet")
		for iind in range(countIndices, countIndices+nInd):
			isUp, nDummy = isDummy(expression[2][iind])
		if nDummy == -1:
			continue
		if expression[1][i] == identity[0]:
			firstList.append(i)
			dummyIndexInfos1.append((isUp, nDummy))
		if expression[1][i] == identity[1]:
			secondList.append(i)
			dummyIndexInfos2.append((isUp, nDummy))
		countIndices += nInd
	if not same:
		firstList,  dummyIndexInfos1 = refineIndices(firstList,  dummyIndexInfos1)
		secondList, dummyIndexInfos2 = refineIndices(secondList, dummyIndexInfos2)
		if not len(firstList) == len(secondList):
			raise Exception("found both objects different number of times")
		for ii in dummyIndexInfos1:
			found = False
			for jj in dummyIndexInfos2:
				if jj[1] == ii[1]:
					if ii[0] == jj[0]:
						raise Exception("two times same index position!!")
					found = True
					break
			if not found:
				raise Exception("unpaired index found!!! OH NO!!!")
		nIdent = len(firstList)
	else:
		checkList, checkInf = refineIndices(firstList,  dummyIndexInfos1)
		if len(checkList) > 0:
			raise Exception("uncontracted indices")
		if not len(firstList)%2 == 0:
			raise Exception("two times same object, but odd number found")
		nIdent = len(firstList)/2
	countIndices = 0
	for i in range(len(expression[1])):
		nInd = expression[1][i].rank
		if countIndices in firstList:
			countIndices += nInd
			continue
		if countIndices in secondList and not same:
			countIndices += nInd
			continue
		indexList = []
		for iind in range(countIndices, countIndices+nInd):
			indexList.append(expression[2][iind])
		countIndices += nInd
		
		retVal *= expression[1][i](*indexList)
	for _ in range(nIdent):
		retVal *= identity[2]
	return retVal


def contractIdentity(identity, expression):
	arg = expression.args
	if  type(arg[-1]).__name__ == "Tuple":
		return applyIndetityToSingeExpression(identity, arg)
	else:
		first = True
		for i in range(len(arg)):
			contracted = contractIdentity(identity, arg[i])
			if first:
				retVal = contracted
				first = False
			else:
				retVal += contracted
		return retVal

if __name__ == "__main__":
	i = getNewLorentzIndex()
	j = getNewLorentzIndex()
	k = getNewLorentzIndex()
	l = getNewLorentzIndex()
	m = getNewLorentzIndex()
	n = getNewLorentzIndex()

	m2 = symbols('m2')

	p = vType('p')
	q = vType('q')

	prood = m2*p(i)*q(j) + 3*m2*p(j)*q(i) + p(k)*p(-k)*q(i)*q(j) + q(k)*q(-k)*p(i)*p(j)
	prod  = m2*p(i)*q(j)
#	print prood
	iden = (p,p,m2)
	print prod
	print contractIdentity(iden, prod)
	print '------------------------------'
	print prood
	print contractIdentity(iden, prood)
	raise Exception
	for su in prood.args:
		print su.args
		for d in su.args[1]:
			print d.name
		for i in su.args[2]:
			print type(i)
			print i.name, i.is_up
