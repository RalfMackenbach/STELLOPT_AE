import  f90nml
import  re
import  time
import  AE_routines         as      ae
import  numpy               as      np
import  matplotlib.pyplot   as      plt
from    scipy.signal        import  savgol_filter
import  multiprocessing     as mp
import warnings
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


################################ DO EXAMPLE CASE ###############################
# If one expands AE in a weakly driven regime, one finds a square dependence on
# the gradients. In a strongly driven regime, one expects a linear dependence.

path = 'gist_files/'
omn_res=100
filter_fraction = 0.1
eta = 0.0
omn_vals = np.logspace(-2,1,omn_res)




filter_window = int(int(omn_res*filter_fraction/2)*2 + 1)

files = sorted(['gist_D3D.txt','gist_HSX.txt','gist_STELLOPT.txt','gist_W7XHM.txt','gist_W7XSC.txt','gist_hsx_norm_00000.dat','gist_hsx_norm_01183.dat'])


if __name__ == "__main__":

    pool = mp.Pool(mp.cpu_count())
    print('Number of cores used: {}'.format(mp.cpu_count()))
    start_time = time.time()
    for file in files:
        if file=='gist_D3D.txt':
            omnigenous = True
        else:
            omnigenous = False
        path_to_file = path+file
        gist_class = gist(path_to_file)
        # time the full integral
        AE_list = pool.starmap(compute_ae_gist, [(gist_class,val,eta*val,omnigenous) for idx, val in np.ndenumerate(omn_vals)])

        AE_arr = np.asarray(AE_list)



        fig, axs = plt.subplots(1, 3)
        c = 0.4
        fig.set_size_inches(c*30, c*12)
        axs[0].loglog( omn_vals, AE_arr )
        log_omns = np.log(omn_vals)
        log_AE   = np.log(AE_arr)
        first_deriv  = list(savgol_filter(log_AE, filter_window, 3,deriv=1)/(log_omns[1]-log_omns[0]))
        second_deriv = list(savgol_filter(log_AE, filter_window, 3,deriv=2)/(log_omns[1]-log_omns[0])**2.0)
        min_value = min(second_deriv)
        min_index = second_deriv.index(min_value)
        crit_grad = omn_vals[min_index]
        print("critical gradient for " + file + " = " + str(omn_vals[min_index]))
        axs[1].plot( omn_vals, first_deriv )
        axs[1].set_xscale("log")
        axs[2].plot( omn_vals, second_deriv )
        axs[2].set_xscale("log")
        axs[2].axvline(x = crit_grad, color = 'r', linestyle='dashed',label ='omn_crit = {}'.format("{:.2e}".format(crit_grad)))
        axs[0].set_xlabel('omn')
        axs[1].set_xlabel('omn')
        axs[2].set_xlabel('omn')
        axs[0].set_ylabel('available energy')
        axs[1].set_ylabel('first derivative')
        axs[2].set_ylabel('second derivative')
        axs[2].legend()
        fig.suptitle(file)
        plt.tight_layout()
        plt.savefig('critgrad_' + file + '.png')
    print("data generated in       --- %s seconds ---" % (time.time() - start_time))
