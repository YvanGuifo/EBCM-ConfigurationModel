import matplotlib.pyplot as plt
import csv
import numpy as np
import EoN
import networkx as nx
import pandas as pd
import random

plt.figure(figsize=(7,7))


step = 0.1     #step size
step_max_plot = 150
xx = np.arange(0, step_max_plot*step, step)

df1 = pd.read_csv('/Users/admin/Desktop/Implementation/MassActionI.csv')
df1 = df1.loc[:, '{#status->#I}']
inc1 = df1[:step_max_plot]

df2 = pd.read_csv('/Users/admin/Desktop/Implementation/CMPowerLawI.csv')
df2 = df2.loc[:, '{#status->#I}']
inc2 = df2[:step_max_plot]

df3 = pd.read_csv('/Users/admin/Desktop/Implementation/CMBimodalI.csv')
df3 = df3.loc[:, '{#status->#I}']
inc3 = df3[:step_max_plot]

df4 = pd.read_csv('/Users/admin/Desktop/Implementation/CMPoissonI.csv')
df4 = df4.loc[:, '{#status->#I}']
inc4 = df4[:step_max_plot]

df5 = pd.read_csv('/Users/admin/Desktop/Implementation/CMHomogeneousI.csv')
df5 = df5.loc[:, '{#status->#I}']
inc5 = df5[:step_max_plot]

plt.plot(xx, inc4, '--', label = 'Poisson', color = 'red')
plt.plot(xx, inc3, '--', label = 'Bimodal', color = 'green')
plt.plot(xx, inc2, '--', label = 'Truncated Power Law', color = 'orange')
plt.plot(xx, inc5, '--', label = 'Homogeneous', color = 'blue')
plt.plot(xx, inc1, '--', label = 'Mass Action', color = 'black')

plt.axis(xmin=0, ymin=0, xmax = 15, ymax = 0.6)
plt.xlabel('$t$')
plt.ylabel('infectious')
#plt.title('Configuration model')
plt.legend()
plt.show()











































