from makeTensors import vType, tType, getNewLorentzIndex, getX, contract, fullContract, getPolarizatrionTensorRest
import sympy
from sympy.utilities.codegen import codegen



if __name__ == "__main__":

	J = 0
	L = 2
	j = 2
	M = 0

	ii = getNewLorentzIndex()
	jj = getNewLorentzIndex()

	mPi = sympy.symbols("mPi")
	p1x, p1y, p1z = sympy.symbols("p1x p1y p1z")
	p2x, p2y, p2z = sympy.symbols("p2x p2y p2z")

	p3x = - p1x - p2x
	p3y = - p1y - p2y
	p3z = - p1z - p2z


	E1 = sympy.sqrt(mPi*mPi + p1x*p1x + p1y*p1y + p1z*p1z)
	E2 = sympy.sqrt(mPi*mPi + p2x*p2x + p2y*p2y + p2z*p2z)
	E3 = sympy.sqrt(mPi*mPi + p3x*p3x + p3y*p3y + p3z*p3z)

	m3pi = E1 + E2 + E3

	P = vType('P')
	P.data = [m3pi, 0, 0, 0]

	p1 = vType("p1")
	p1.data = [E1, p1x, p1y, p1z]

	p2 = vType("p2")
	p2.data = [E2, p2x, p2y, p2z]

	p3 = vType("p3")
	p3.data = [E3, p3x, p3y, p3z]

	GP = tType("GP")
	GP.data = [	[0, 0, 0, 0],
			[0,-1, 0, 0],
			[0, 0,-1, 0],
			[0, 0, 0,-1]]

	g = tType("g")
	g.data = [	[1, 0, 0, 0],
			[0,-1, 0, 0],
			[0, 0,-1, 0],
			[0, 0, 0,-1]]

	p12 = p1(ii) + p2(ii)
	m122 = (p12(ii)*p12(-ii)).data

	gp = g(ii,jj) - p12(ii)*p12(jj)/m122

	K = (p1(ii) + p2(ii) - p3(ii))/2
	KP = GP(ii,jj)*K(-jj)

	k = (p1(ii) - p2(ii))/2
	kp = gp(ii,jj)*k(-jj)

	XL = getX(KP, GP, L)
	Xj = getX(kp, gp, j)
	
	T = contract(XL, Xj, J, P)
	eps = getPolarizatrionTensorRest(J,M)
	
	try:
		ampl = 	fullContract(T,eps).data
	except AttributeError:
		ampl = T*eps
	print("simpl")
	ampl = sympy.simplify(ampl)
	print(ampl)
	




