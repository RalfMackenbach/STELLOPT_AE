import  f90nml
import  AE_routines     as  ae
import  re
import  numpy as np
import  matplotlib.pyplot   as plt
from    matplotlib          import cm
import  itertools
import  tikzplotlib
from matplotlib import rc
import matplotlib.lines as mlines
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import  matplotlib          as mpl


class gist:
    def __init__(self, path_to_file):
        self.nml = f90nml.read(path_to_file)
        # Import entire txt file and assign to self
        txt_RGX       = (open(path_to_file, 'r')).read()
        self.txt_file = txt_RGX

        # Import all columns
        pattern_functions  = re.compile(r"""(?<=\/)(?s)(.*$)""",flags=re.MULTILINE|re.DOTALL)
        x = (re.findall(pattern_functions, self.txt_file)[0])
        x_split = x.split('\n ')

        l = []
        for item in x_split:
            subl = []
            for num in item.split():
                subl.append(float(num))
            l.append(subl)
        l_new = list2 = [x for x in l if x != []]
        self.functions = np.asarray(l_new)




def compute_ae_gist(gist_class,omn,omt):
    # import relevant parameters from gist class
    params = gist_class.nml['parameters']
    q0 = params['q0']
    n_pol = params['n_pol']
    gridpoints = params['gridpoints']

    # make variious arrays needed for calculations
    theta_arr   = np.linspace(0.0, n_pol*2*np.pi, gridpoints)
    B           = gist_class.functions[:,3]
    jac         = gist_class.functions[:,4]
    dBdx        = gist_class.functions[:,5]
    dBdy        = gist_class.functions[:,6]
    sqrt_g      = 1/jac
    Delta_x     = q0
    Delta_y     = q0
    B_ave       = np.trapz(q0*B*B*sqrt_g,theta_arr)/np.trapz(q0*B*sqrt_g,theta_arr)
    b_arr       = B/B_ave
    dbdx_arr    = dBdx/B_ave
    dbdy_arr    = dBdy/B_ave

    # set scalars
    dlnTdx = -omt
    dlnndx = -omn
    L_tot  = np.trapz(q0*b_arr*sqrt_g,theta_arr)

    # set numerical parameters
    z_min = 1e-4
    z_max = 40
    z_res = 10000
    lam_res = 10000
    Delta_theta = 1e-10
    del_sing = 1e02
    ae_val = ae.ae_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr,lam_res,z_res,z_min,z_max,Delta_theta,del_sing,L_tot)
    return ae_val


############################### MAKE PRL PLOT #################################
# Note that this plot is slightly different than the one in the PRL, since we
# are using slightly different numerical methods here. For example, in the PRL
# code I explicitly calculated arclength(theta), and then proceeded to
# execute various integrals as f d(arclength). Here, however, I execute the
# integrals as f d(arclength)/dtheta dtheta. The factor d(arclength)/dtheta can
# be written out in terms of the jacobian and magnetic field strength. Both
# codes should behave identically in the limit of very fine theta grids.




path = 'gist_files/'



Q_list=np.asarray([7.645847154751877,37.11184344194595,66.28563923274653,88.89918001896879,0.2214852174780842,0.8211983903213222,1.848364581960216,8.272401253718314,0.42191750695843805,0.31121588405080763,0.9969773188039313,5.042615883328997,0.9204760920288282,79.72394526587243,2272.082920266666])
name_list=["D3D - omt=0.0, omn=1.0","D3D - omt=0.0, omn=2.0","D3D - omt=0.0, omn=3.0","D3D - omt=0.0, omn=4.0","HSX - omt=0.0, omn=1.0","HSX - omt=0.0, omn=2.0","HSX - omt=0.0, omn=3.0","HSX - omt=0.0, omn=4.0","W7XHM - omt=0.0, omn=1.0","W7XHM - omt=0.0, omn=2.0","W7XHM - omt=0.0, omn=3.0","W7XHM - omt=0.0, omn=4.0","W7XSC - omt=0.0, omn=3.0","W7XSC-Tdrive - omt=3.0, omn=0.0","W7XSC-Tdrive - omt=3.0, omn=0.0"]
print(name_list)

AE_list = []

for omn in range(1,5):
    file = 'gist_D3D.txt'
    path_to_file = path+file
    gist_class = gist(path_to_file)
    AE_val = compute_ae_gist(gist_class,omn,0.0)
    AE_list = np.append(AE_list,AE_val)

for omn in range(1,5):
    file = 'gist_HSX.txt'
    path_to_file = path+file
    gist_class = gist(path_to_file)
    AE_val = compute_ae_gist(gist_class,omn,0.0)
    AE_list = np.append(AE_list,AE_val)

for omn in range(1,5):
    file = 'gist_W7XHM.txt'
    path_to_file = path+file
    gist_class = gist(path_to_file)
    AE_val = compute_ae_gist(gist_class,omn,0.0)
    AE_list = np.append(AE_list,AE_val)

file = 'gist_W7XSC.txt'
path_to_file = path+file
gist_class = gist(path_to_file)
AE_val = compute_ae_gist(gist_class,3.0,0.0)
AE_list = np.append(AE_list,AE_val)

file = 'gist_W7XSC.txt'
path_to_file = path+file
gist_class = gist(path_to_file)
AE_val = compute_ae_gist(gist_class,0.0,3.0)
AE_list = np.append(AE_list,AE_val)

file = 'gist_D3D.txt'
path_to_file = path+file
gist_class = gist(path_to_file)
AE_val = compute_ae_gist(gist_class,0.0,3.0)
AE_list = np.append(AE_list,AE_val)

AE_list = np.asarray(AE_list)


# calculate and print power laws
TT_list = np.asarray([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 7])
z,cov = np.polyfit(np.log(AE_list), np.log(Q_list/(TT_list**(5/2))), 1, cov=True)
print("exponent and SD for Q = AE^n")
print(z[0],np.sqrt(cov[0,0]))
print("lin term and SD for Q = AE^n")
print(z[1],np.sqrt(cov[1,1]))


marker = itertools.cycle(('P', 'P', 'P', 'P',
                          'o', 'o', 'o', 'o',
                          '^', '^', '^', '^',
                          'v','^','P'))
cmap = cm.get_cmap('plasma_r', 5)
all_vals = cmap([0,1,2,3,4])
colormap = ListedColormap(all_vals[1:5])

color  = itertools.cycle((1, 2, 3, 4,
                          1, 2, 3, 4,
                          1, 2, 3, 4,
                          3))





plt.style.use("seaborn-bright")
fig = plt.figure(figsize=(3.375, 2.75))
ax  = fig.gca()
ax.set_yscale('log')
ax.set_xscale('log')
# ax.grid(which='minor', alpha=0.1,zorder=0)
ax.grid(which='minor', color='0.95', zorder=0)
ax.grid(which='major', color='0.90', zorder=0)
z,cov = np.polyfit(np.log(AE_list), np.log(Q_list/(TT_list**(5/2))), 1, cov=True)
polynomial = np.poly1d(z)
log_y_fit = polynomial(np.log(AE_list))
plt.plot(AE_list, np.exp(log_y_fit), '-',color='black',zorder=1,linewidth=1.0)

for idx, AE_val in enumerate(AE_list):
    if (TT_list[idx]==1.0):
        color_val = ((next(color)-1)/3 - 1.5)*1.0 + 1.5
        color_plt = colormap(color_val)
    else:
        color_plt='grey'
    ax.scatter(float(AE_list[idx]), float(Q_list[idx])/float((TT_list[idx]**(5/2))), marker=next(marker), color=color_plt, zorder=2, s=20)


DIIID       = mlines.Line2D([], [], color='black', marker='P', linestyle='None',
                          markersize=5.5, label='DIII-D')
HSX         = mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                          markersize=5.5, label='HSX')
W7XHM       = mlines.Line2D([], [], color='black', marker='^', linestyle='None',
                          markersize=5.5, label='W7-X (HM)')
W7XSC       = mlines.Line2D([], [], color='black', marker='v', linestyle='None',
                          markersize=5.5, label='W7-X (SC)')
NCSX        = mlines.Line2D([], [], color='black', marker='s', linestyle='None',
                          markersize=5.5, label='NCSX')
ax.legend(handles=[DIIID, HSX, W7XHM, W7XSC])
norm = mpl.colors.Normalize(vmin=.5,vmax=4.5)
cbar = fig.colorbar(cm.ScalarMappable(norm=norm,cmap=colormap),ax=ax,ticks=[1, 2, 3, 4],label="density gradient")
cbar.ax.tick_params(size=0)
ax.xaxis.set_tick_params(which='major', direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', direction='in', top='on')
ax.yaxis.set_tick_params(which='major', direction='in', top='on')
ax.yaxis.set_tick_params(which='minor', direction='in', top='on')
plt.xlabel("available energy")
plt.ylabel("turbulent energy flux")
ax.set_axisbelow(True)
fig.tight_layout()
plt.subplots_adjust(left=0.15, right=1.0, top=0.99, bottom=0.13)
plt.show()
