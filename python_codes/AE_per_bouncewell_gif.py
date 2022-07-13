import  f90nml
import  re
import  time
import  numpy               as      np
import  AE_routines         as      ae
import  matplotlib.pyplot   as      plt
import  matplotlib.colors   as      mplc
import  multiprocessing     as      mp
from    matplotlib          import  cm
import  numpy.ma            as      ma
import  matplotlib          as      mpl
from    matplotlib          import  rc
from    scipy.interpolate   import  interp1d
import  warnings
import  imageio
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
# rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
## for Palatino and other serif fonts use:
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
warnings.filterwarnings("ignore")

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





def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,**kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw))
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.2, -0.09),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=270,**kw))
                     for text,color in zip(list_of_strings[::-1],list_of_colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(1.10, 0.3),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)



def plot_AE_per_bouncewell(theta_arr,b_arr,dbdx_arr,lam_arr,bw,ae_list,ae,omn,idx,ae_arr,omn_arr):
    c = 0.5
    # shift by pi
    fig ,ax = plt.subplots(1,2)
    fig.set_size_inches(c*30, c*12)
    ax[1].set_xlim(min(theta_arr),max(theta_arr))

    list_flat = []
    for val in ae_list:
        list_flat.extend(val)

    max_val = max(list_flat)
    cm_scale = lambda x: x
    colors_plot = [cm.plasma(cm_scale(np.asarray(x) * 1.0/max_val)) for x in ae_list]

    # iterate over all values of lambda
    for idx_lam, lam in enumerate(lam_arr):
        b_val = 1/lam

        # iterate over all bounce wells
        for idx_bw, _ in enumerate(ae_list[idx_lam]):
            # check if well crosses boundary
            if(bw[idx_lam][idx_bw][0] > bw[idx_lam][idx_bw][1]):
                ax[1].plot([bw[idx_lam][idx_bw][0], max(theta_arr)], [b_val, b_val], color=colors_plot[idx_lam][idx_bw])
                ax[1].plot([min(theta_arr), bw[idx_lam][idx_bw][1]], [b_val, b_val], color=colors_plot[idx_lam][idx_bw])
            # if not normal plot
            else:
                ax[1].plot([bw[idx_lam][idx_bw][0], bw[idx_lam][idx_bw][1]], [b_val, b_val], color=colors_plot[idx_lam][idx_bw])

    ax[1].plot(theta_arr,b_arr,color='black',linewidth=2)
    ax2 = ax[1].twinx()
    ax2.plot(theta_arr, dbdx_arr, 'red')
    ax2.plot(theta_arr,theta_arr*0.0,linestyle='dashed',color='red')
    ax2.spines['right'].set_color('red')
    ax2.tick_params(axis='y', colors='red',labelcolor='red',direction='in')
    ax[1].set_xlabel(r'$\theta$')
    ax[1].set_xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
    ax[1].set_xticklabels([r'$-\pi$', r'$-\pi/2$',r'$0$', r'$\pi/2$', r'$\pi$'])
    ax[1].tick_params(axis='both',direction='in')
    plt.colorbar(cm.ScalarMappable(norm=mplc.Normalize(vmin=0.0, vmax=max_val/ae, clip=False), cmap=cm.plasma), ax=ax[1],location='bottom',label=r'$A_\lambda/A$') #'%.3f'
    multicolor_ylabel(ax2,[r'$\partial B / \partial x $', r'$B$,'],['black','red'],'y',0.1)
    ax[0].loglog(omn_arr,ae_arr)
    ax[0].axvline(x = omn, color = 'black', linestyle='dashed')
    ax[0].set_xlabel(r'$\omega_n$')
    ax[0].set_ylabel(r'$A$')
    plt.savefig('ae_per_bw_{}.png'.format("{:03d}".format(idx)),format='png',dpi=200)
    plt.close()
    return 'ae_per_bw_{}.png'.format("{:03d}".format(idx))


def compute_ae_gist(gist_class,omn,omt,omnigenous):
    # import relevant parameters from gist class
    params = gist_class.nml['parameters']
    q0 = np.float64(params['q0'])
    n_pol = params['n_pol']
    gridpoints = params['gridpoints']

    # make various arrays needed for calculations
    theta_arr   = np.float64(np.linspace(0.0, n_pol*2*np.pi, gridpoints))
    B           = np.float64(gist_class.functions[:,3])
    jac         = np.float64(gist_class.functions[:,4])
    dBdx        = np.float64(gist_class.functions[:,5])
    dBdy        = np.float64(gist_class.functions[:,6])
    sqrt_g      = 1/jac
    B_ave       = np.trapz(q0*B*B*sqrt_g,theta_arr)/np.trapz(q0*B*sqrt_g,theta_arr)
    b_arr       = B/B_ave
    dbdx_arr    = dBdx/B_ave
    dbdy_arr    = dBdy/B_ave
    Delta_x     = q0
    Delta_y     = q0

    # set scalars
    dlnTdx = np.float64(-omt)
    dlnndx = np.float64(-omn)
    L_tot  = np.trapz(q0*b_arr*sqrt_g,theta_arr)

    # set numerical parameters
    lam_res = 1000
    Delta_theta = 1e-10
    del_sing = 0.00
    ae_val = ae.ae_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr,lam_res,Delta_theta,del_sing,L_tot,omnigenous)
    return ae_val



def compute_ae_over_z_gist(gist_class,omn,omt,omnigenous,idx,ae_arr,omn_arr):
    # import relevant parameters from gist class
    params = gist_class.nml['parameters']
    q0 = params['q0']
    n_pol = params['n_pol']
    gridpoints = params['gridpoints']

    # make variious arrays needed for calculations
    theta_arr   = np.linspace(0.0, n_pol*2*np.pi, gridpoints) - np.pi
    B           = gist_class.functions[:,3]
    jac         = gist_class.functions[:,4]
    dBdx        = gist_class.functions[:,5]
    dBdy        = gist_class.functions[:,6]
    sqrt_g      = 1/jac
    B_ave       = np.trapz(q0*B*B*sqrt_g,theta_arr)/np.trapz(q0*B*sqrt_g,theta_arr)
    b_arr       = B/B_ave
    dbdx_arr    = dBdx/B_ave
    dbdy_arr    = dBdy/B_ave
    Delta_x     = q0
    Delta_y     = q0

    # set scalars
    dlnTdx = -omt
    dlnndx = -omn
    L_tot  = np.trapz(q0*b_arr*sqrt_g,theta_arr)

    # set numerical parameters
    lam_res = 500
    Delta_theta = 1e-10
    del_sing = 0.0
    bw, lam_arr, ae_list, ae_tot = ae.ae_total_over_z(q0,dlnTdx,dlnndx,Delta_x,Delta_y,b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr,lam_res,Delta_theta,del_sing,L_tot,omnigenous)
    b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr = ae.make_per(b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr,Delta_theta)
    filename = plot_AE_per_bouncewell(theta_arr,b_arr,dbdx_arr,lam_arr,bw,ae_list,ae_tot,omn,idx,ae_arr,omn_arr)
    return filename






################################ DO EXAMPLE CASE ###############################
# Do AE per trapping well

path = 'gist_files/'
file = 'gist_D3D.txt'
omn_res=400
omn_vals = np.logspace(-2,1,omn_res)
eta = 0.0


if file=='gist_D3D.txt':
    omnigenous = True
else:
    omnigenous = False

gist_class = gist(path + file)
filenames = []
images = []



if __name__ == "__main__":

    pool = mp.Pool(mp.cpu_count())
    print('Nu24mber of cores used: {}'.format(mp.cpu_count()))
    start_time = time.time()
    AE_list = pool.starmap(compute_ae_gist, [(gist_class,val,eta*val,omnigenous) for idx, val in np.ndenumerate(omn_vals)])
    AE_list = np.asarray(AE_list)
    filename = pool.starmap(compute_ae_over_z_gist, [(gist_class,omn,eta*omn,omnigenous,idx[0],AE_list,omn_vals)  for idx, omn in np.ndenumerate(omn_vals)])
    filenames.extend(filename)
    print("data generated in       --- %s seconds ---" % (time.time() - start_time))
    pool.close()

    print("generating gif...")
    for filename in sorted(filenames):
        print("adding " + str(filename) + "...")
        images.append(imageio.imread(filename))
    imageio.mimsave(file+'.mp4', images, fps=24)
    print("done!")
