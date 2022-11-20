import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from sympy import symbols, Eq, solve, exp

K = 1
r0 = 0.05
betar = 1 
time = 6720
v0=1
betadelta=1
delta0 =.02
zeta0 = .03
zeta02 = .01
betazeta = 1
betazeta2 = 5
amax = 2
betaa = 0.5 
vd = 2
ft=0
k = .01
ther_start = 0
ther_end = time
txchange = [ther_start]
assess = [ther_start]
min_summ = []
txdose = 1

IC = [0.01,1]

length = 168
vbound = 1.2 #.8,1.2, 1.6; 1.9
tmbound = 0.7 #.3, .5, 7
curr = 0


def grow(v):
    return r0*(v/v0)**betar

def death(v):
    return delta0*math.exp(betadelta*(v-v0))
    
def chemo(v):
    return zeta0*math.exp(-betazeta*(v-v0))

def evoLV(X,t):    
            
    C = X[0]
    v = X[1]
    tm = v*C
    
    def ther(): #total marker
        global curr
        if t>=ther_start and t<=ther_end:
            if t>assess[-1]+length:
                assess.append(t)
                if tm>tmbound:
                    if curr!=1:
                        txchange.append(t)
                    curr = 1
                else:
                    if curr != 0:
                        txchange.append(t)
                    curr = 0
            return curr
        else:
            return 0
        
    
    # def ther(): #mean marker
    #     global curr
    #     if t>=ther_start and t<=ther_end:
    #         if t>assess[-1]+length:
    #             assess.append(t)
    #             if v>vbound:
    #                 if curr!=1:
    #                     txchange.append(t)
    #                 curr = 1
    #             else:
    #                 if curr != 0:
    #                     txchange.append(t)
    #                 curr = 0
    #         return curr
    #     else:
    #          return 0
        

    def G(v):
        return grow(v)*(1-C/K)-death(v)-txdose*ther()*chemo(v)
    
    dCdt = C*G(v)
            
    dvdt = k*((betar*r0/v0*(v/v0)**(betar-1))*(1-C/K)-grow(v)*C/K*(-betaa*(amax-1)/amax)-(betadelta*death(v))-txdose*ther()*(-betazeta*chemo(v)))
    
    dxvdt = np.array([dCdt, dvdt])

    return dxvdt

intxv = np.array(IC)
time_sp = np.linspace(0,time,time*100-1)
pop = odeint(evoLV, intxv,time_sp,atol = 1e-13, rtol=1e-13)
print((pop[:,1][-1]))

    
plt.figure()
plt.subplot(311)
plt.title('Chemotherapy: Daily Visits, Low Threshold')
plt.plot(time_sp,pop[:,0],lw=3,color='b')
#plt.axhline(y=.35, color='c',lw=1, linestyle='dashed')
plt.xticks([])
ax = plt.gca()
txchange.append(ther_end)
for i in range(len(txchange)-1):
    if i%2==1:
        ax.axvspan(txchange[i], txchange[i+1], facecolor='grey')

plt.xlim(0,time)
plt.ylabel('Pop Size')

plt.subplot(312)
ax = plt.gca()
for i in range(len(txchange)-1):
    if i%2==1:
        ax.axvspan(txchange[i], txchange[i+1], facecolor='grey')
plt.plot(time_sp,pop[:,1],lw=3,color='b')
plt.xticks([])
plt.ylabel('Mean Marker')
plt.xlim(0,time)

plt.subplot(313)
ax = plt.gca()
for i in range(len(txchange)-1):
    if i%2==1:
        ax.axvspan(txchange[i], txchange[i+1], facecolor='grey')
plt.plot(time_sp,pop[:,0]*pop[:,1],lw=3,color='b')
plt.axhline(y=tmbound, color='c',lw=1, linestyle='dashed')
plt.ylabel('Total Marker')
plt.xlim(0,time)

plt.xticks([0,1344,4032,6720],[0,8,24,40])
plt.xlabel('Weeks')
plt.show()
