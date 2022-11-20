import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from math import exp

K = 1
r0 = 0.05
betar = 1 #1,2
v0 = 1
betadelta = 1
delta0 = .02
zeta0 = .03*4 #*4
zeta02 = .01*4 #*4
betazeta = 1
betazeta2 = 5
amax = 2
betaa = 0.5 #0.5,0
vd = 2
ft = 0
k = .01
chemo_start = 1344 #8 weeks
chemo_end = chemo_start+2688 #24 weeks
inter_start = chemo_start+2688
inter_end = inter_start+2688
time=inter_end
min_summ = []

total_dose = 160
num_doses = 16 
txdose = 1

IC = [.6,1]

def grow(v):
    return r0*(v/v0)**betar

def death(v):
    return delta0*math.exp(betadelta*(v-v0))

def chemo(v):
    return zeta0*math.exp(-betazeta*(v-v0))

def target(v):
    return zeta02*((v/v0)**betazeta2)

def fc(t):
    if t>=chemo_start and t<=chemo_end:
        txtime = total_dose/num_doses/txdose
        period = (chemo_end-chemo_start)/num_doses
        if (t-chemo_start)%period < txtime:
            return txdose
        else:
            return 0
    else:
        return 0
    
def ft(t):
    if t>=inter_start and t<=inter_end:
        txtime = total_dose/num_doses/txdose
        period = (inter_end-inter_start)/num_doses
        if (t-inter_start)%period < txtime:
            return txdose
        else:
            return 0
    else:
        return 0

def evoLV(X, t):    
    
    C = X[0]
    v = X[1]
    
    def G(v):
        return grow(v)*(1-C/K)-death(v)-fc(t)*chemo(v)-ft(t)*target(v)
    
    dCdt = C*G(v)
        
    dvdt = k*((betar*r0/v0*(v/v0)**(betar-1))*(1-C/K)-grow(v)*C/K*(-betaa*(amax-1)/amax)-(betadelta*death(v))-fc(t)*(-betazeta*chemo(v))-ft(t)*(betazeta2*zeta02/v0*(v/v0)**(betazeta2-1)))
        
    dxvdt = np.array([dCdt, dvdt])

    return dxvdt

intxv = np.array(IC)
time_sp = np.linspace(0,time,num=time*100-1)
pop = odeint(evoLV, intxv,time_sp,atol=1e-12,rtol=1e-12)
plt.figure()
plt.subplot(211)
plt.title('Sequential Therapy: High Tx Intensity, Low Evolv')
plt.plot(time_sp,pop[:,0],lw=3,color='b')
plt.xlim(0,inter_end)
plt.xticks([])
ax = plt.gca()
ax.axvspan(chemo_start, chemo_end, facecolor=(0.5,0.5,0.5))
ax.axvspan(inter_start, inter_end, facecolor='r') #firebrick
plt.ylabel('Pop Size, x')
plt.subplot(212)
ax = plt.gca()
ax.axvspan(chemo_start, chemo_end, facecolor=(0.5,0.5,0.5))
ax.axvspan(inter_start, inter_end, facecolor='r')
plt.plot(time_sp,pop[:,1],lw=3,color='b')
plt.ylabel('Biomarker Level, v')
plt.xlabel('Weeks')
plt.xticks([0,1344,4032,6720],[0,8,24,40])
plt.xlim(0,inter_end)
plt.show()
print(min(pop[:,0]))
