import numpy as n
import pylab as plab
import scipy.integrate as integrate

#Choroid Plexus AA
def dap(ap,ac):
    return p*ap*(1-ap/kappa)-c*ap*(1-ac/lagr)-dea*ap

#Blood Brain Barrier DHA
def ddp(dp,dc):
    return b*dp*(1-dp/xi)-e*dp*(1-dc/chi)-dd*dp

#CSF AA
def dac(ac,ap,aa,an):
    return (1-ac/lagr)*(c*ap+s*aa)-alpha*f*ac*(1-an/upsilon)-dea*ac

#CSF DHA
def ddc(dc,dp,dn,da):
    return (1-dc/chi)*(e*dp+g*dn+z*da)-h*dc*(1-dn/omega)-i*dc*(1-da/psi)-dd*dc

#Neuronal AA
def dan(an,ac,dn):
    return (1-an/upsilon)*(alpha*f*ac+j*dn)-alpha*k*an*(1-dn/omega)-dea*an

#Neuronal DHA
def ddn(dn,dc,an):
    return (1-dn/omega)*(h*dc+alpha*k*an)-j*dn*(1-an/upsilon)-g*dn*(1-dc/chi)-dd*dn

#Neuronal ROS
def drn(rn,dn,an):
    return alpha*q-alpha*r*k*an*rn*(1-dn/omega)-dr*rn

#Astrocyte AA
def daa(aa,da,ac):
    return (1-aa/tau)*m*da-s*aa*(1-ac/lagr)-dea*aa

#Astrocyte DHA
def dda(da,dc,aa):
    return (1-da/psi)*i*dc-m*da*(1-aa/tau)-z*da*(1-dc/chi)-dd*da

def derivs(state, t):

    ap,dp,ac,dc,an,dn,rn,aa,da = state
    deltaap = dap(ap,ac)
    deltadp = ddp(dp,dc)
    deltaac = dac(ac,ap,aa,an)
    deltadc = ddc(dc,dp,dn,da)
    deltaan = dan(an,ac,dn)
    deltadn = ddn(dn,dc,an)
    deltarn = drn(rn,dn,an)
    deltaaa = daa(aa,da,ac)
    deltada = dda(da,dc,aa)
    return deltaap, deltadp, deltaac, deltadc, deltaan, deltadn, deltarn, deltaaa, deltada

#Values for model constants: obtained from original ODE system

dea=0.017 #Overall brain loss of ascorbic acid is >2%
dd=0.017 #Assuming a similar rate for degradation of DHAA
dr = 0.02 #Death of ROS due to SOD
q=5 #Change to make appropriate, relatively arbitrary
r=1 #1 to 1 relation between ROS production and AA_N to DHAA_N
p = .15 #CP calculation: fit to Michaelis-Menten kinetics uptake model (neglected diffusion b/c it's impact was negligible) [Lam, Daniel 1985]
kappa = 60 #NIH data
b = 0.25 #V_blood brain barrier estimation
c=0.72 #CP calculation
e=0.72 #V estimation
f=0.1 #CP calculation
g=0.72 #V estimation
h=0.72 #V estimation
i=0.1 #GLUT scaling for lack of GLUT3 receptors [Maher 1994]
j=0.9 #"This is very rapid" [Pinto paper]
k=.8 #Estimate myself
m=0.9 #"Very rapid"
s=0.1 #CP calculation
z=0.1 #GLUT scaling
lagr=160 #NIH data - The following are all just the normal, physiological parameters.
xi=20 #NIH data
chi=30 #Estimte myself
tau=3000 #Estimate myself
alpha=.2 #Change to test different physiological conditions, arbitrary again
psi=200 #NIH data
upsilon=7000 #NIH data, estimate myself
omega=100 #Estimate myself


#Initial Concentrations
ap0=10
dp0=5
ac0=12
dc0=8
an0=20
dn0=5
rn0=10
aa0=15
da0=7

t = n.arange(0, 200, 0.1)

y0 = [ap0,dp0,ac0,dc0,an0,dn0,rn0,aa0,da0]  #The initial [intra, plasma] state vector
y = integrate.odeint(derivs, y0, t)
ap = y[:, 0]
dp = y[:, 1]
ac = y[:, 2]
dc = y[:, 3]
an = y[:, 4]
dn = y[:, 5]
rn = y[:, 6]
aa = y[:, 7]
da = y[:, 8]

plab.figure()
plab.plot(t, ap, label='Asc Plasma')
plab.plot(t, dp, label='DHA Blood Brain')
plab.plot(t, ac, label='Asc CSF')
plab.plot(t, dc, label='DHA CSF')
plab.plot(t, an, label='Asc Neuron')
plab.plot(t, dn, label='DHA Neuron')
plab.plot(t, rn, label='ROS')
plab.plot(t, aa, label='Asc Astrocyte')
plab.plot(t, da, label='DHA Astrocyte')
plab.xlabel('Time')
plab.ylabel('Concentration')
plab.title('Vitamin C Trajectories')
plab.grid()
plab.legend()

plab.figure()
plab.plot(an,rn)
plab.xlabel('Asc Neuron')
plab.ylabel('ROS')
plab.title('Phase Plane')

#Make direction field plot with quiver
# apmax = 1.1 * ap.max()
# dpmax = 1.1 * dp.max()
# acmax = 1.1 * ac.max()
# dcmax = 1.1 * dc.max()
# anmax = 1.1 * an.max()
# dnmax = 1.1 * dn.max()
# rnmax = 1.1 * rn.max()
# aamax = 1.1 * aa.max()
# damax = 1.1 * da.max()
#
# AP,DP,AC,DC,AN,DN,RN,AA,DA = n.meshgrid(n.arange(-1, apmax), n.arange(-1, dpmax), n.arange(-1, acmax), n.arange(-1, dcmax), n.arange(-1, anmax), n.arange(-1, dnmax), n.arange(-1, rnmax), n.arange(-1, aamax), n.arange(-1, damax))
# dAP = dap(AP,AC)
# dDP = ddp(DP,DC)
# dAC = dac(AC,AP,AA,AN)
# dDC = ddc(DC,DP,DN,DA)
# dAN = dan(AN,AC,DN)
# dDN = ddn(DN,DC,AN)
# dRN = drn(RN,AN)
# dAA = daa(AA,DA,AC)
# dDA = dda(DA,DC,AA)
# plab.quiver(AP,DP,AC,DC,AN,DN,RN,AA,DA, dAP, dDP, dAC, dDC, dAN, dDN, dRN, dAA, dDA)
#
#
# AP,DP,AC,DC,AN,DN,RN,AA,DA = n.meshgrid(n.arange(-1, apmax,.1), n.arange(-1, dpmax,.1), n.arange(-1, acmax,.1), n.arange(-1, dcmax,.1), n.arange(-1, anmax,.1), n.arange(-1, dnmax,.1), n.arange(-1, rnmax,.1), n.arange(-1, aamax,.1), n.arange(-1, damax,.1))
# dAP = dap(AP,AC)
# dDP = ddp(DP,DC)
# dAC = dac(AC,AP,AA,AN)
# dDC = ddc(DC,DP,DN,DA)
# dAN = dan(AN,AC,DN)
# dDN = ddn(DN,DC,AN)
# dRN = drn(RN,AN)
# dAA = daa(AA,DA,AC)
# dDA = dda(DA,DC,AA)
#
# plab.contour(AP, AC, dAP, levels=[0], linewidths=3, colors='black')
# plab.contour(DP,DC, dDP, levels=[0], linewidths=3, colors='black')
# plab.contour(AC,AP,AA,AN,dAC, levels=[0], linewidths=3, colors='black')
# plab.contour(DC,DP,DN,DA,dDC, levels=[0], linewidths=3, colors='black')
# plab.contour(AN,AC,DN,dAN, levels=[0], linewidths=3, colors='black')
# plab.contour(DN,DC,AN,dDN, levels=[0], linewidths=3, colors='black')
# plab.contour(RN,AN,dRN, levels=[0], linewidths=3, colors='black')
# plab.contour(AA,DA,AC,dAA, levels=[0], linewidths=3, colors='black')
# plab.contour(DA,DC,AA,dDA, levels=[0], linewidths=3, colors='black')
# plab.title('Trajectory, Direction Field, and Nullclines')

plab.show()