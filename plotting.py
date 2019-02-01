from examples import *
import matplotlib.pyplot as plt
from numpy import linspace, array, amin
from mpl_toolkits.mplot3d import Axes3D

#curve = get_hilbert_curve()
#curve = get_peano_curve()
#curve = get_meurthe_curve()
#curve = get_coil_curve()
#curve = get_serpentine_curve()
#curve = get_R_curve()
#curve = get_haverkort_curve_1()
curve = get_haverkort_curve_2()
#curve = get_tokarev_curve()
#curve = get_discontinuous_curve()
#curve = get_morton_curve()

def plot_curve(subdiv_n, dim, genus, numb_sub):
    
    subdiv_n = array(subdiv_n)
    subdiv_n = subdiv_n/(genus**((numb_sub+1)/dim))
    
    for k in range(dim):
        subdiv_n[:,k] = subdiv_n[:,k] - amin(subdiv_n[:,k]) + 1/(2*genus**((numb_sub+1)/(dim)))

    if dim == 2:

        ticks = linspace(0,1,genus**(numb_sub/2)+1)
    
        plt.figure()
        plt.gcf().set_size_inches(9,9)
        plt.xticks(ticks,[])
        plt.yticks(ticks,[])
        plt.grid(True)
        plt.axis([0,1,0,1])
        plt.plot(subdiv_n[:,0],subdiv_n[:,1],'k')

    elif dim == 3:
    
        ticks = linspace(0,1,genus**(numb_sub/3)+1)
    
        fig = plt.figure(figsize = (9,9))
        ax = fig.gca(projection='3d')
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_zlim(0,1)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)
        ax.tick_params(colors = 'w')
        ax.plot(subdiv_n[:,0],subdiv_n[:,1],subdiv_n[:,2],'k')

numb_sub = 1
plot_curve(curve.get_subdivision(k=numb_sub).proto, curve.dim, curve.div**curve.dim, numb_sub)
