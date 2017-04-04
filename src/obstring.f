      subroutine obstring(obs,string)
c     this subroutine returns a string describing the function of obs
c     
c     author  : Steffen A. Bass
c     date    : 10.08.93
c     revision: 0.9
c
      IMPLICIT REAL*8(A-H,O-Z)

      integer obs
      character*15 string
c
c
      if(obs.gt.100)goto 2000
c now goto the desired observable
c           1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
      goto( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,99,12,13,14,15,16,17,18,19,20,
     &     21,22,23,24,25,26,27,28,29,30,99,99,99,99,35,36,37,38,39,40,
     &     41,42,99,99,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &     61,62,63,64,65,66,67,68,99,99,99,99,99,99,99,99,99,99,99,99,
     &     99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,301
     & ), obs
c
c now goto event-observables
 2000 continue
c          101 102 103 104 105 106 107 108 109 110 111 112 113 114 115
      goto(  1,102,103,104,105,106,107, 99, 99, 99,111,112,113,114,115,
     &     116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,
     &     131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,
     &     146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,
     &     161,162,163,164,165,166,167,168,169,170, 99, 99, 99, 99, 99,
     &      99, 99, 99, 99, 99, 99,182,183,184,185,186,187, 20,189, 99,
     &      99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99
     & ), (obs-100)
c here come the observables
c feel free to add your favourite observable;
c in case of changes, please contact the authors 
c for a general revision update.
c
 1    string='dN/d(x/yobs)'
      return
 2    string='p_x [GeV]'
      return
 3    string='p_y [GeV]'
      return
 4    string='p_z [GeV]'
      return
 5    string='p_t [GeV]'
      return
 6    string='Y'
      return
 7    string='E_tot [GeV]'
      return
 8    string='phi [deg]'
      return
 9    string='theta [deg]'
      return
 10   string='|p| [GeV]'
      return
 12   string='mass m [GeV]'
      return
 13   string='m_t [GeV]'
      return
 14   string='BMS: R(I,5)'
      return
 15   string='3*el.charge'
      return
 16   string='PDG_ID'
      return
 17   string='# hard coll'
      return
 18   string='# space bra'
      return
 19   string='# time bra'
      return
 20   string=' v2'
      return
 21   string=' virt. [GeV^2]'
      return
 22   string='time [fm/c]'
      return
 23   string='x [fm]'
      return
 24   string='y [fm]'
      return
 25   string='z [fm]'
      return
 26   string='Y/Y_proj'
      return
 27   string='p_x/p_proj'
      return
 28   string='p_t/p_proj'
      return
 29   string='x_Bj'
      return
 30   string='qzz'
      return
 35   string='p_x/A [GeV]'
      return
 36   string='p_y/A [GeV]'
      return
 37   string='p_z/A [GeV]'
      return
 38   string='p_t/A [GeV]'
      return
 39   string='m_t-m_0 [GeV]'
      return
 40   string='part. index'
      return
 41   string='E_T'
      return 
 42   string='trans. press.'
      return
 45   string='t_freezeout [fm/c] '
      return
 46   string='r_x(fr) [fm]'
      return
 47   string='r_y(fr) [fm]'
      return
 48   string='r_z(fr) [fm]'
      return
 49   string='|r(fr)| [fm]'
      return
 50   string='|r_t(fr)| [fm]' 
      return
 51   string='|r_t| [fm]' 
      return
 52   string='r_t*p_t/|p_t|'
      return
 53   string='eta '
      return
 54   string='log(x_Bj)'
      return
 55   string='Y_init'
      return
 56   string='(p_z)^2'
      return
 57   string='(p_t)^2'
      return
 58   string='medium T'
      return
 59   string='initial p_t'
      return
 60   string='initial E'
      return
 61   string='initial particle weight'
      return
 62   string='v_4'
      return
 63   string='v4/v2^2'
      return
 64   string='vT(cell)'
      return
 65   string='v2(cell)'
      return
 66   string='vx(cell)'
      return
 67   string='vy(cell)'
      return
 68   string='vz(cell)'
      return
 69   string='vx^2-vy^2(cell)'
      return
 301  string='sp-tree'
      return
 99   string='jodeldidodeldoe'
      return
 102  string='b [fm]'
      return
 103  string='part Mult(Grp)'
      return
 104  string='beam ener. (lab)'
      return
 105  string='sqrt(s)'
      return
 106  string='beam mom. (lab)'
      return
 107  string='beam mom. (CM)'
      return
 111  string='(beta_ref-cm)x'
      return
 112  string='(beta_ref-cm)y'
      return
 113  string='(beta_ref-cm)z'
      return
 114  string='<px> CM-sys'
      return
 115  string='<py> CM-sys'
      return
 116  string='<pz> CM-sys'
      return
 117  string='sigma(px)'
      return
 118  string='sigma(py)'
      return
 119  string='sigma(pz)'
      return
 120  string='EV1 mom-flow'
      return
 121  string='theta_flow_1'
      return
 122  string='phi_EV1_RP'
      return
 123  string='EV2 mom-flow'
      return
 124  string='theta_flow_2'
      return
 125  string='phi_EV2_RP'
      return
 126  string='EV3 mom-flow'
      return
 127  string='theta_flow_3'
      return
 128  string='phi_EV3_RP'
      return
 129  string='rat:EV1/EV3'
      return
 130  string='rat:EV2/EV3'
      return
 131  string='EV1 sphericity'
      return
 132  string='theta_sph_1'
      return
 133  string='phi_EV1_RP_sph'
      return
 134  string='EV2 sphericity'
      return
 135  string='theta_sph_2'
      return
 136  string='phi_EV2_RP_sph'
      return
 137  string='EV3 sphericity'
      return
 138  string='theta_sph_3'
      return
 139  string='phi_EV3_RP_sph'
      return
 140  string='rat:EV1/EV3_s'
      return
 141  string='rat:EV2/EV3_s'
      return
 142  string='EV1 kin-flow'
      return
 143  string='theta_kin_1'
      return
 144  string='phi_EV1_RP_kin'
      return
 145  string='EV2 kin-flow'
      return
 146  string='theta_kin_2'
      return
 147  string='phi_EV2_RP_kin'
      return
 148  string='EV3 kin-flow'
      return
 149  string='theta_kin_3'
      return
 150  string='phi_EV3_RP_kin'
      return
 151  string='rat:EV1/EV3_k'
      return
 152  string='rat:EV2/EV3_k'
      return
 153  string='<cos(theta_f)>'
      return
 154  string='<cos(theta_s)>'
      return
 155  string='<cos(theta_k)>'
      return
 156  string='gamma_ref_CM'
      return
 157  string='<px,dir>'
      return
 158  string='sigma(px,dir)'
      return
 159  string='<|p|> CM-sys'
      return
 160  string='sigma(|p|)'
      return
 161  string='<E_kin> CMS'
      return
 162  string='sigma(E_kin)'
      return
 163  string='<p_t> CMS'
      return
 164  string='sigma(p_t)'
      return
 165  string='p_XYT CMS'
      return
 166  string='sig(beta_x)'
      return
 167  string='sig(beta_y)'
      return
 168  string='sig(beta_z)'
      return
 169  string='QZZ/p^2 CMS'
      return
 170  string='# part(flw)'
      return
 182  string='# of coll'
      return
 183  string='<rt>'
      return
 184  string='<tf>'
      return
 185  string='<pt>'
      return
 186  string='E_t'
      return
 187  string='E_tot in ev.'
      return
 189  string='Ecc_x'
      return
      end












