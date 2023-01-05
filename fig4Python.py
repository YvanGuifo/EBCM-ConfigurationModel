import networkx as nx
import EoN
from collections import defaultdict
import matplotlib.pyplot as plt
import scipy
import numpy as np
import random
import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd

plt.figure(figsize=(8,7))

colors = ['#5AB3E6','#FF2000','#009A80','#E69A00', '#CD9AB3', '#0073B3','#F0E442', '#7f00ff', '#0000ff']
rho = 0.025
Nbig=500000
Nsmall = 5000
tau =0.6
gamma = 1.

def poisson():
    return np.random.poisson(5)
def PsiPoisson(x):
    return np.exp(-5*(1-x))
def DPsiPoisson(x):
    return 5*np.exp(-5*(1-x))

bimodalPk = {8:0.5, 2:0.5}
def PsiBimodal(x):
    return (x**8 +x**2)/2.
def DPsiBimodal(x):
    return(8*x**7 + 2*x)/2.

def homogeneous():
    return 5
def PsiHomogeneous(x):
    return x**5
def DPsiHomogeneous(x):
    return 5*x**4


PlPk = {}
exponent = 1.418184432
kave = 0
for k in range(1,81):
    PlPk[k]=k**(-exponent)*np.exp(-k*1./40)
    kave += k*PlPk[k]

normfact= sum(PlPk.values())
for k in PlPk:
    PlPk[k] /= normfact

def PsiPowLaw(x):
    #print PlPk
    rval = 0
    for k in PlPk:
        rval += PlPk[k]*x**k
    return rval

def DPsiPowLaw(x):
    rval = 0
    for k in PlPk:
        rval += k*PlPk[k]*x**(k-1)
    return rval


def get_G(N, Pk):
    while True:
        ks = []
        for ctr in range(N):
            r = random.random()
            for k in Pk:
                if r<Pk[k]:
                    break
                else:
                    r-= Pk[k]
            ks.append(k)
        if sum(ks)%2==0:
            break
    G = nx.configuration_model(ks)
    return G

def process_degree_distribution(Gbig, color, Psi, DPsi,  label):

    N= Gbig.order()#N is arbitrary, but included because our implementation of EBCM assumes N is given.
    t, S, I, R = EoN.EBCM(N, lambda x: (1-rho)*Psi(x), lambda x: (1-rho)*DPsi(x), tau, gamma, 1-rho)
    plt.plot(t, I/N, color = color, label=label)
    
#Erdos Renyi
Gbig = nx.fast_gnp_random_graph(Nbig, 5./(Nbig-1))
process_degree_distribution(Gbig, colors[1], PsiPoisson, DPsiPoisson, r'Poisson')

#Bimodal
Gbig = get_G(Nbig, bimodalPk)
process_degree_distribution(Gbig, colors[2], PsiBimodal, DPsiBimodal, r'Bimodal')


#Powerlaw
Gbig = get_G(Nbig, PlPk)
process_degree_distribution(Gbig, colors[3], PsiPowLaw, DPsiPowLaw, r'Truncated Power Law')

#Homogeneous
Gbig = get_G(Nbig, {5:1.})
process_degree_distribution(Gbig, colors[8], PsiHomogeneous, DPsiHomogeneous, r'Homogeneous')

step = 0.1     #step size
step_max_plot = 150

xx = np.arange(0, step_max_plot*step, step)

df1 = pd.read_csv('/usr/admin/Desktop/Implementation/MassActionI.csv')
df1 = df1.loc[:, '{#status->#I}']
inc1 = df1[:step_max_plot]

plt.plot(xx, inc1, label = 'Mass Action', color = 'black')

plt.axis(xmin=0, ymin=0, xmax = 15, ymax = 0.6)
plt.xlabel('$t$')
plt.ylabel('infectious')
plt.legend()
plt.show()
