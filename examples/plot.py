# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 16:08:29 2021

@author: rred_
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 18:10:50 2016

@author: hezhang
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:16:54 2016

@author: hezhang
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:34:13 2016

@author: hezhang
"""

import numpy as np
import matplotlib.pyplot as plt 

# set global settings
def init_plotting():
    plt.rcParams['figure.figsize'] = (9,8)
    plt.rcParams['font.size'] = 24
    plt.rcParams['font.family'] = 'Times New Roman'
#    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
#    plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
    plt.rcParams['lines.markersize'] = 10
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.loc'] = 'center left'
    plt.rcParams['legend.numpoints'] = 1
    plt.rcParams['axes.linewidth'] = 1
#    plt.gca().spines['right'].set_color('none')
#    plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')


filename_trunc = 'track_truncated_map_cosy'
filename_symp = 'track_symplectic_map_type_1_cosy'

# filename_sim = 'runjspec2/output_imp_1kv_400ns-01-2'
# filename_exp = 'Max/jspec_test/rmsbunchlength_1kv_400ns'
data_trunc = np.loadtxt(filename_trunc+'.txt', dtype='float64', skiprows=0)
x_t = data_trunc[:,1]
xp_t = data_trunc[:,2]
y_t = data_trunc[:,3]
yp_t = data_trunc[:,4]

data_symp = np.loadtxt(filename_symp+'.txt', dtype='float64', skiprows=0)
x_s = data_symp[:,1]
xp_s = data_symp[:,2]
y_s = data_symp[:,3]
yp_s = data_symp[:,4]


#
fig, ax1 = plt.subplots()
#ax1.set_xlim(0,60)
# ax1.set_ylim(0, 2.2)
#
n_every = 1
line_1 = ax1.plot(x_t[::n_every],xp_t[::n_every],'o', color = 'g', label=r'truncated map', linewidth=1)
line_2 = ax1.plot(x_s[::n_every],xp_s[::n_every],'.', color = 'r', label=r'symplectic map', linewidth=1)
ax1.set_xlabel(r'x')
ax1.set_ylabel(r'xp' )

###
plt.grid()
lines = line_1 + line_2 

labels = [l.get_label() for l in lines]
plt.legend(lines, labels, loc="upper left")
plt.savefig("x-xp.pdf", format="pdf")
plt.show()

fig, ax1 = plt.subplots()
#ax1.set_xlim(0,60)
# ax1.set_ylim(0, 2.2)
#
n_every = 1
line_1 = ax1.plot(y_t[::n_every],yp_t[::n_every],'o', color = 'g', label=r'truncated map', linewidth=1)
line_2 = ax1.plot(y_s[::n_every],yp_s[::n_every],'.', color = 'r', label=r'symplectic map', linewidth=1)
ax1.set_xlabel(r'y')
ax1.set_ylabel(r'yp' )

###
plt.grid()
lines = line_1 + line_2 

labels = [l.get_label() for l in lines]
plt.legend(lines, labels, loc="upper right")
plt.savefig("y-yp.pdf", format="pdf")
plt.show()
