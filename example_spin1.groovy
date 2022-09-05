//Mikhail Mikhasenko 
//1.02.2016

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

  def gp = 'gp_ab = g_ab-P_a*P_b/(P_c*P^c)'.t
  def kp = gp >> 'kp_a = k^b*gp_ab'.t  

  def m1rel = 'k1_a*k1^a = m**2'.t, m2rel = 'k2_a*k2^a = m**2'.t, m3rel = 'k3_a*k3^a = m**2'.t
  def s12rel = '2*k1_a*k2^a = s12-2*m**2'.t
  def s23rel = '2*k2_a*k3^a = s23-2*m**2'.t
  def s31rel = '2*k3_a*k1^a = s31-2*m**2'.t
  def delta = 'd_a^a = 4'.t
  def invariants = m1rel & m2rel & m3rel & s12rel & s23rel & s31rel & delta

  def X2gen = orbT('_a'.si)
  def scal2 = 'P_c*P^c = s12'.t & 'P_a = k1_a+k2_a'.t & 'k_a = 1/2*(k1_a-k2_a)'.t
  X2 = (kp & gp & ExpandAndEliminate[scal2] & ExpandAndEliminate) >> X2gen
  X2 <<= (ExpandAndEliminate[invariants] & ExpandAndEliminate)
  println X2.toString(LaTeX)
  println "\$\\\\[1cm]\$"

  def scal1 = 'P_c*P^c = s'.t & 'P_a = k1_a+k2_a+k3_a'.t & 'k_a = 1/2*(k1_a+k2_a-k3_a)'.t
  X1 = (kp & gp & ExpandAndEliminate[scal1] & ExpandAndEliminate) >> X2gen
  X1 <<= (ExpandAndEliminate[invariants] & ExpandAndEliminate)
  println X1.toString(LaTeX)
  println "\$\\\\[1cm]\$"

  def M = ('d^a_y'.t * X1) * ('g^ay'.t * X2)
  M <<= (ExpandAndEliminate & 
	 ExpandAndEliminate[invariants & 's31 = s+3*m**2-s12-s23'.t] & 
	 ExpandAndEliminate)
  println M.toString(LaTeX)


}
