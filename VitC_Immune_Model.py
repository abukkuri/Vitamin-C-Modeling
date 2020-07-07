import numpy as n
import pylab as plab
import scipy.integrate as integrate

#Vitamin C
def dC(C,V):
    return c0-dc*C-a*V/VM

#Phacogyctes
def dP(P,C,V):
    return b*P*(1-P/PM)-P*(dp+(e+f*C/CM)*V/VM)

#Host Cells
def dH(H,C,P,V):
    return g*H*(1-H/HM)-H*(dh+h*(i-j*C/CM)*P/PM+k*V/VM)

#Viruses
def dV(V,C,P):
    return (l-m*C/CM)*V*(1-V/VM)-V*(dv+(z+p*C/CM)*(P/PM)*(q+(r+s*C/CM)+(u+v*C/CM)))

def derivs(state, t):

    C,P,H,V = state
    deltaC = dC(C,V)
    deltaP = dP(P,C,V)
    deltaH = dH(H,C,P,V)
    deltaV = dV(V,C,P)
    return deltaC, deltaP, deltaH, deltaV

#Values for model constants: obtained from original ODE system

c0 = 0 #50
dc = .05
a = .1
b = 1.5
PM = 400
dp = .1
e = .5
f = 0.1
g = 1.5
HM = 1000
dh = .1
h = .2
i = .6
j = .3
k = .4
l = 2.5
m = 1.25
VM = 300
dv = .1
z = .6
p = .3
q = .2
r = .6
s = .2
u = .2
v = .2
CM = 2000

#Initial Concentrations
C0 = 2000
P0 = 300
H0 = 1000
V0 = 500

t = n.arange(0, 200, 0.1)

y0 = [C0,P0,H0,V0]
y = integrate.odeint(derivs, y0, t)
C = y[:, 0]
P = y[:, 1]
H = y[:, 2]
V = y[:, 3]

plab.figure()
plab.plot(t, C, label='Vitamin C')
plab.plot(t, P, label='Phagocytic Cells')
plab.plot(t, H, label='Host Cells')
plab.plot(t, V, label='Virus')
plab.xlabel('Time')
plab.ylabel('Quantity')
plab.title('Immune System Dynamics under Infection')
plab.grid()
plab.legend()

plab.figure()
plab.plot(C,V)
plab.xlabel('Vitamin C')
plab.ylabel('Virus')
plab.title('Phase Plane')

plab.show()