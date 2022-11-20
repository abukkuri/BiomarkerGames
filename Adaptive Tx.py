import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from sympy import symbols, Eq, solve, exp

K = 1
r0 = 0.05
betar = 1 #1,2
time = 6720
v0=1
betadelta=1
delta0 =.02
zeta0 = .03 #*4
zeta02 = .01 #*4
betazeta = 1
betazeta2 = 5
amax = 2
betaa = 0.5 #.5,0
k = .01 #.1
ther_start = 1344
ther_end = time
txchange = [ther_start]
assess = [ther_start]
min_summ = []
txdose = 1

IC = [0.6,1]

length = 168
vbound = 1.9 #.8,1.2, 1.6; 1.9
tmbound = 0.68 #.3, .45, .6; .65
curr = [1,0]


def grow(v):
    return r0*(v/v0)**betar

def death(v):
    return delta0*math.exp(betadelta*(v-v0))
    
def chemo(v):
    return zeta0*math.exp(-betazeta*(v-v0))

def target(v):
    return zeta02*((v/v0)**betazeta2)

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
                    if curr!=[0,1]:
                        txchange.append(t)
                    curr = [0,1]
                else:
                    if curr != [1,0]:
                        txchange.append(t)
                    curr = [1,0]
            return curr
        else:
            return [0,0]
        
    
    # def ther(): #mean marker
    #     global curr
    #     if t>=ther_start and t<=ther_end:
    #         if t>assess[-1]+length:
    #             assess.append(t)
    #             if v>vbound:
    #                 if curr!=[0,1]:
    #                     txchange.append(t)
    #                 curr = [0,1]
    #             else:
    #                 if curr != [1,0]:
    #                     txchange.append(t)
    #                 curr = [1,0]
    #         return curr
    #     else:
    #         return [0,0]
        

    def G(v):
        return grow(v)*(1-C/K)-death(v)-txdose*ther()[0]*chemo(v)-txdose*ther()[1]*target(v)
    
    dCdt = C*G(v)
            
    dvdt = k*((betar*r0/v0*(v/v0)**(betar-1))*(1-C/K)-grow(v)*C/K*(-betaa*(amax-1)/amax)-(betadelta*death(v))-txdose*ther()[0]*(-betazeta*chemo(v))-txdose*ther()[1]*(betazeta2*zeta02/v0*(v/v0)**(betazeta2-1)))
    
    dxvdt = np.array([dCdt, dvdt])

    return dxvdt

intxv = np.array(IC)
time_sp = np.linspace(0,time,time*100-1)
pop = odeint(evoLV, intxv,time_sp,atol = 1e-12, rtol=1e-12)
print(min(pop[:,0]))


# #Threshold Total Marker
# arr1 = [0.17607357549527136,0.1477651461243536,0.11693061543594616,0.10051454332766453,0.05581283627524448,0.031319954288481444,0.012794858677358512,0.005395540393704198,0.0008348716038133836,0.00041963023050223334,0.00022588191830465634,0.0001368764718977938,8.117878650086665e-05,3.666751596154172e-05,1.8860183741170112e-05,6.577492587736533e-06,0.17272368823946713]
# times1 = [.2,.25,.3,.35,.4,.45,.5,.55,.6,.61,.62,.63,.64,.65,.66,.67,.68]

# #Threshold Mean Marker
# arr = [0.17607357549527136,0.022270030119410426,0.032116421996958876,0.04592154505244598,0.06271997051340132,0.08648775710761475,0.09931156517810734,0.11842203900340605,0.13667822034107507,0.14506854854836884,0.15818365331970202,0.1477651461243536,0.11693061543594616,0.10051454332766453,0.06021839330733347,0.04889436046367521,0.02684088657480923,0.012794858677358512,0.005395540393704198,0.0027095127401553965,0.0008348716038133836,0.00032558876576835186,8.117878650086665e-05,2.5647662365596796e-05,5.265706870310936e-06,0.17272368823946713]
# times = [.7,.75,.8,.85,.9,.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95]


# plt.figure()
# plt.title('Effect of Marker Thresholds on Adaptive Therapy')
# plt.ylabel('Nadir of Cancer Population')
# plt.xlabel('Marker Thresholds')
# plt.plot(times,arr,lw=3,color='b')
# plt.plot(times1,arr1,lw=3,color='r')
# ax = plt.gca()
# ax.set_yscale('log')

    
plt.figure()
plt.subplot(211)
plt.title('Adaptive Therapy Mean Marker: High Evolvability')
plt.plot(time_sp,pop[:,0],lw=3,color='b')
#plt.axhline(y=.35, color='c',lw=1, linestyle='dashed')
plt.xticks([])
ax = plt.gca()
txchange.append(ther_end)
for i in range(len(txchange)-1):
    if i%2==0:
        ax.axvspan(txchange[i], txchange[i+1], facecolor='grey')
    else:
        ax.axvspan(txchange[i], txchange[i+1], facecolor='r')

plt.xlim(0,time)
plt.ylabel('Pop Size, x')
plt.subplot(212)
ax = plt.gca()
    
for i in range(len(txchange)-1):
    if i%2==0:
        ax.axvspan(txchange[i], txchange[i+1], facecolor='grey')
    else:
        ax.axvspan(txchange[i], txchange[i+1], facecolor='r')
plt.plot(time_sp,pop[:,1],lw=3,color='b')
plt.axhline(y=1.19749, color='c',lw=1, linestyle='dashed')
plt.ylabel('Drug Resistance, v')
plt.xlim(0,time)
plt.xticks([0,1344,4032,6720],[0,8,24,40])
plt.xlabel('Weeks')
plt.show()

# #Tmbound, vbound, Length, Tx Dose

# '''for i in np.arange(.5,2,.05):
#     vbound = i
#     pop = odeint(evoLV, intxv,time_sp,atol = 1e-12, rtol=1e-12)
#     a = min(pop[:,0])
#     if a < 0:
#         a = 0
#     min_summ.append([i,a])
    
# min_summ = np.array(min_summ)
# plt.figure()
# plt.title('Effect of Threshold on Uninformed Therapy')
# plt.ylabel('Nadir of Cancer Population')
# plt.xlabel('Tx Switching Threshold')
# plt.plot(min_summ[:,0],min_summ[:,1],lw=3,color='k')'''