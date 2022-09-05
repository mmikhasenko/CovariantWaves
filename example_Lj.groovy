//Mikhail Mikhasenko
//16.02.2016

@Grab(group = 'cc.redberry', module = 'groovy', version = '1.1.10')
import cc.redberry.groovy.Redberry

import cc.redberry.core.indices.*
import cc.redberry.core.indexgenerator.IndexGeneratorImpl

import static cc.redberry.groovy.RedberryPhysics.*
import static cc.redberry.groovy.RedberryStatic.*
import static cc.redberry.core.tensor.Tensors.*
import static cc.redberry.core.context.OutputFormat.*

use(Redberry){

  setSymmetric 'gp_ab'
  setAntiSymmetric 'eps_abcd'
  def orbT
  orbT = {ll -> 
    //if(!(t instanceof Indices)) return
    int dim = ll.size()
    //construct tensors reqursively
    if( dim == 0) return '1'.t
    if( dim == 1) return ('kp'+ll[0].toStringIndex()).t
    if( dim >= 2) {
      def fpart1 = '0'.t
      for(int i=0; i<=dim-1; i++) {
	def Xi = orbT(ll[(0..dim-1).findAll{ it!=i }])
	fpart1 += Xi * (('kp'+ll[i].toStringIndex()).t)
      }
      //second part
      IndexGeneratorImpl generator = new IndexGeneratorImpl(ll.toArray());
      def fpart2 = '0'.t
      for(int i=0; i<dim-1; i++) {
	for(int j=i+1; j<dim; j++) {
	  def newind = generator.generate(ll[i].type);
	  def ind = ll[(0..dim-1).findAll{ (it!=i)&&(it!=j) }]
	  ind = ind+newind
	  def Xij = orbT(ind)
	  fpart2 += Xij * (('gp'+ll[i].toStringIndex()+ll[j].toStringIndex()).t) * (('kp'+newind.invert().toStringIndex()).t)
	}
      }

      def fpart = ('(2*'+dim+'-1)/'+dim+'**2').t * fpart1 - ('2/'+dim+'**2').t * fpart2
      fpart = ExpandAndEliminate >> fpart
      //fpart = 'gp_ab = gp_ba'.t >> fpart
      fpart = 'kp^a*gp_ab = kp_b'.t >> fpart
      fpart = 'kp^a*gp_ca = kp_c'.t >> fpart
      return fpart
    }
  }


  def contract = {T1, T2, j ->
    def iT1 = T1.indices.free, iT2 = T2.indices.free;
    int i1 = iT1.size(),       i2 = iT2.size()
    if(j>i1+i2 || j<(i1-i2)) return
    def fexpr = T1*T2;
    if((i1+i2-j)%2==0) {
      //contruct
      int k = (i1+i2-j)/2
      if(k==0) return fexpr;
      if(k!=0) {
	for(int r=0;r<=k-1;r++) {
	  fexpr *= ('g'+iT1[r].invert().toStringIndex()+iT2[r].invert().toStringIndex()).t
	}
      }
      return fexpr;
    } else if((i1+i2-j)%2==1) {
      //contruct
      int k = (i1+i2-j-1)/2
      if(k!=0) {
	for(int r=0;r<=k-1;r++) {
	  fexpr *= ('g'+iT1[r].invert().toStringIndex()+iT2[r].invert().toStringIndex()).t
	}
      }
      IndexGeneratorImpl generator = new IndexGeneratorImpl(iT1+iT2);      
      //generator.s
      def newind1 = generator.generate(iT1[0].type)
      def newind2 = generator.generate(iT2[0].type)
      def fstr = 'eps'+iT1[k].toUpper().toStringIndex()+iT2[k].toUpper().toStringIndex() +
	newind1.toUpper().toStringIndex() + newind2.toStringIndex() +
	' * P'+newind1.toStringIndex()
	//println fstr;
      fexpr *= ( fstr ).t
      return fexpr
    }
    return T1*T2
  }

//  def test = contract(orbT('_abc'.si),orbT('_rty'.si),1)
//  test <<= ExpandAndEliminate
//  println test

  //println orbT('_ab'.si)


  int L = 3
  int j = 1
  int J = 2

  //orbital AngMom
  def indL = []; for(int i=0;i<L;i++) indL<<i
  def XL = orbT(indL as Indices)

  //j of isobar
  def indj = []; for(int i=0;i<j;i++) indj<<i+L
  def Xj = orbT(indj as Indices)

  //************************************************************//
  //************************************************************//
  //************************************************************//

  //ortogonality
  def gp = 'gp_ab = g_ab-P_a*P_b/(P_c*P^c)'.t
  def kp = gp >> 'kp_a = k^b*gp_ab'.t  
  //invariants
  def m1rel = 'k1_a*k1^a = m**2'.t, m2rel = 'k2_a*k2^a = m**2'.t, m3rel = 'k3_a*k3^a = m**2'.t
  def s12rel = '2*k1_a*k2^a = s12-2*m**2'.t
  def s23rel = '2*k2_a*k3^a = s23-2*m**2'.t
  def s31rel = '2*k3_a*k1^a = s31-2*m**2'.t
  def delta = 'd_a^a = 4'.t
  def invariants = m1rel & m2rel & m3rel & s12rel & s23rel & s31rel & delta

  def scalj = 'P_c*P^c = s12'.t & 'P_a = k1_a+k2_a'.t & 'k_a = 1/2*(k1_a-k2_a)'.t
  Xj <<= (kp & gp & ExpandAndEliminate[scalj] & ExpandAndEliminate)
  Xj <<= (ExpandAndEliminate[invariants] & ExpandAndEliminate)
  println Xj.toString(LaTeX)
  println "\$\\\\[1cm]\$"

  def scalL = 'P_c*P^c = s'.t & 'P_a = k1_a+k2_a+k3_a'.t & 'k_a = 1/2*(k1_a+k2_a-k3_a)'.t
  XL <<= (kp & gp & ExpandAndEliminate[scalL] & ExpandAndEliminate)
  XL <<= (ExpandAndEliminate[invariants] & ExpandAndEliminate)
  println XL.toString(LaTeX)
  println "\$\\\\[1cm]\$"

  //************************************************************//
  //************************************************************//
  def M = contract(Xj,XL,J)
  M <<= 'P_a = k1_a+k2_a+k3_a'.t
  M <<= (ExpandAndEliminate & 
	 ExpandAndEliminate[invariants & 's31 = s+3*m**2-s12-s23'.t] & 
	 ExpandAndEliminate)
  println M.toString(LaTeX)

}
