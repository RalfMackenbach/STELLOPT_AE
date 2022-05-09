import  f90nml
import  re
import  time
import  numpy               as      np
import  AE_routines         as      ae
import  matplotlib.pyplot   as      plt
import  matplotlib.colors   as      mplc
import  multiprocessing     as      mp
from    matplotlib          import  cm
import  warnings
import  imageio
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



def plot_AE_per_bouncewell(theta_arr,b_arr,dbdx_arr,lam_arr,bw,ae_list,ae,omn,idx):

    fig ,ax = plt.subplots()
    ax.set_xlim(min(theta_arr),max(theta_arr))

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
                ax.plot([bw[idx_lam][idx_bw][0], max(theta_arr)], [b_val, b_val], color=colors_plot[idx_lam][idx_bw])
                ax.plot([min(theta_arr), bw[idx_lam][idx_bw][1]], [b_val, b_val], color=colors_plot[idx_lam][idx_bw])
            # if not normal plot
            else:
                ax.plot([bw[idx_lam][idx_bw][0], bw[idx_lam][idx_bw][1]], [b_val, b_val], color=colors_plot[idx_lam][idx_bw])

    ax.plot(theta_arr,b_arr,color='black',linewidth=2)
    ax2 = ax.twinx()
    ax2.plot(theta_arr,dbdx_arr,color='crimson')
    ax2.plot(theta_arr,theta_arr*0.0,linestyle='dotted',color='crimson')
    ax2.tick_params(axis='y', labelcolor='crimson')
    plt.title('omn={}'.format("{:.2e}".format(omn)))
    plt.colorbar(cm.ScalarMappable(norm=mplc.Normalize(vmin=0.0, vmax=max_val/ae, clip=False), cmap=cm.plasma), ax=ax,location='bottom',label='A_lam/A') #'%.3f'
    plt.savefig('ae_per_bw_{}.png'.format("{:03d}".format(idx)),dpi=200)
    plt.close()
    return 'ae_per_bw_{}.png'.format("{:03d}".format(idx))



def compute_ae_over_z_gist(gist_class,omn,omt,omnigenous,idx):
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
    lam_res = 200
    Delta_theta = 1e-10
    del_sing = 0.0
    bw, lam_arr, ae_list, ae_tot = ae.ae_total_over_z(q0,dlnTdx,dlnndx,Delta_x,Delta_y,b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr,lam_res,Delta_theta,del_sing,L_tot,omnigenous)
    b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr = ae.make_per(b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr,Delta_theta)
    filename = plot_AE_per_bouncewell(theta_arr,b_arr,dbdx_arr,lam_arr,bw,ae_list,ae_tot,omn,idx)
    return filename






################################ DO EXAMPLE CASE ###############################
# Do AE per trapping well

path = 'gist_files/'
file = 'gist_HSX.txt'
omn_res=100
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
    print('Number of cores used: {}'.format(mp.cpu_count()))
    start_time = time.time()
    filename = pool.starmap(compute_ae_over_z_gist, [(gist_class,omn,eta*omn,False,idx[0])  for idx, omn in np.ndenumerate(omn_vals)])
    filenames.extend(filename)
    print("data generated in       --- %s seconds ---" % (time.time() - start_time))
    pool.close()

    print("generating gif...")
    kargs = { 'duration': 0.1 }
    for filename in sorted(filenames):
        print("adding " + filename + "...")
        images.append(imageio.imread(filename))
    imageio.mimsave(file+'.gif', images, **kargs)
    print("done!")
