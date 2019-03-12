from examples import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sub_numb = 2

#curve = get_hilbert_curve()
#curve = get_peano_curve()
#curve = get_meurthe_curve()
#curve = get_coil_curve()
#curve = get_serpentine_curve()
#curve = get_R_curve()
curve = get_haverkort_curve_1()
#curve = get_haverkort_curve_2()
#curve = get_tokarev_curve()
#curve = get_discontinuous_curve()
#curve = get_morton_curve()

def plot_curve(subdiv_n, dim, genus, sub_numb):
    
    #Масштабируем координаты, делим все значения на genus**((sub_numb+1)/dim)
    subdiv_n = [[i/(genus**((sub_numb+1)/dim)) for i in subdiv_n[j]] for j in range(len(subdiv_n))]   
    #Находим минимум по каждой коодинате x,y,z b и т.д., 
    #это необходимо для кривых, у которых точки входа и выхода находятся не в вершинах 
    max_col = [min([row[i] for row in subdiv_n]) for i in range(dim)]
    #Сдвиг всех координат на 1/(2*genus**((sub_numb+1)/(dim))) с учетом минимальных значений по координатам
    subdiv_n = [[subdiv_n[j][i] + max_col[i] + 1/(2*genus**((sub_numb+1)/(dim))) for i in range(dim)] for j in range(len(subdiv_n))]
    
    #Функция создания сетки для 
    def linspace(genus,sub_numb):
        ticks = list(range(0,int(genus**(sub_numb/2)+1)))
        ticks = [i/genus**(sub_numb/2) for i in ticks]
        return ticks
    
    if dim == 2:

        ticks = linspace(genus,sub_numb)
    
        plt.figure()
        plt.gcf().set_size_inches(9,9)
        plt.xticks(ticks,[])
        plt.yticks(ticks,[])
        plt.grid(True)
        plt.axis([0,1,0,1])
        sub1 = [row[0] for row in subdiv_n]
        sub2 = [row[1] for row in subdiv_n]
        plt.plot(sub1,sub2,'k')

    elif dim == 3:
    
        ticks = linspace(genus,sub_numb)
    
        fig = plt.figure(figsize = (9,9))
        ax = fig.gca(projection='3d')
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_zlim(0,1)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)
        ax.tick_params(colors = 'w')
        sub1 = [row[0] for row in subdiv_n]
        sub2 = [row[1] for row in subdiv_n]
        sub3 = [row[2] for row in subdiv_n]
        ax.plot(sub1,sub2,sub3,'k')


plot_curve(curve.get_subdivision(k=sub_numb).proto, curve.dim, curve.div**curve.dim, sub_numb)
