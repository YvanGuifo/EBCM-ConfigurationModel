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
meank = 5
tau =0.6
gamma = 1.

bimodalPk = {8:0.5, 2:0.5}
def PsiBimodal(x):
    return (x**8 +x**2)/2.
def DPsiBimodal(x):
    return(8*x**7 + 2*x)/2.

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

def process_degree_distribution(Gbig, color, Psi, DPsi, symbol, label):
    #G = get_G(N, Pk)
    t, S, I, R = EoN.fast_SIR(Gbig, tau, gamma, rho=rho)
    plt.plot(t, I*1./Gbig.order(), color = color)
    
    N= Gbig.order()#N is arbitrary, but included because our implementation of EBCM assumes N is given.
    t, S, I, R = EoN.EBCM(N, lambda x: (1-rho)*Psi(x), lambda x: (1-rho)*DPsi(x), tau, gamma, 1-rho)
    plt.plot(t, I/N, symbol, color = color, label=label)
    
#Bimodal
Gbig = get_G(Nbig, bimodalPk)
actual_degrees = [d for v, d in Gbig.degree()]
actual_degrees = list(actual_degrees)
square_degrees = [(i**2)-i for i in actual_degrees]
sum_square_degre = sum(square_degrees)
## Mean square degree
mean_square_degre = sum_square_degre/(N*meank)
mean_square_degre
