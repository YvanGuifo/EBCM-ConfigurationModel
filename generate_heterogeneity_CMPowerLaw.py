# Compute (⟨k^2- k⟩/⟨k⟩)

import networkx as nx
import EoN
from collections import defaultdict
import matplotlib.pyplot as plt
import scipy
import numpy as np
import random

colors = ['#5AB3E6','#FF2000','#009A80','#E69A00', '#CD9AB3', '#0073B3','#F0E442']
rho = 0.025
rho1 = 0.001
Nbig=500000
N = 500000
tau =0.6
gamma = 1.
meank = 5

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
    
#Power Law
Gbig = get_G(Nbig, PlPk)
actual_degrees = [d for v, d in Gbig.degree()]
actual_degrees = list(actual_degrees)
square_degrees = [(i**2)-i for i in actual_degrees]
sum_square_degre = sum(square_degrees)
## Mean square degree
mean_square_degre = sum_square_degre/(N*meank)
mean_square_degre
