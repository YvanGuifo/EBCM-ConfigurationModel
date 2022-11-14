import networkx as nx
import EoN
from collections import defaultdict
import matplotlib.pyplot as plt
import scipy
import numpy as np
import random

colors = ['#5AB3E6','#FF2000','#009A80','#E69A00', '#CD9AB3', '#0073B3','#F0E442']
rho = 0.025
meanK = 5
Nbig=500000
N = 500000
tau =0.6
gamma = 1.

def homogeneous():
    return 5
def PsiHomogeneous(x):
    return x**5
def DPsiHomogeneous(x):
    return 5*x**4
  
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
    
#Homogeneous
Gbig = get_G(Nbig, {5:1.})
actual_degrees = [d for v, d in Gbig.degree()]
actual_degrees = list(actual_degrees)
square_degrees = [(i**2)-i for i in actual_degrees]
sum_square_degre = sum(square_degrees)
## Mean square degree
mean_square_degre = sum_square_degre/(N*meanK)
mean_square_degre
