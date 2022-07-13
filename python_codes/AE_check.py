import  f90nml
import  AE_routines         as  ae
import  re
import  numpy as np
import  matplotlib.pyplot   as plt
import  matplotlib          as mpl
from    matplotlib          import  rc
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
# rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
## for Palatino and other serif fonts use:
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

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
    lam_res = 100
    Delta_theta = 1e-10
    del_sing = 0.0
    ae_val, as_val = ae.ae_total_MJG(q0,dlnTdx,dlnndx,Delta_x,Delta_y,b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr,lam_res,Delta_theta,del_sing,L_tot)
    return ae_val, as_val




################################ DO EXAMPLE CASE ###############################
# If one expands AE in a weakly driven regime, one finds a square dependence on
# the gradients. In a strongly driven regime, one expects a linear dependence.
# Check if that is what you find! This picture changes for devices which have
# omega_psi = 0

path = 'gist_files/'
omn_res=100
omn_hlf=int(omn_res/2)
omn_vals = np.logspace(-6,6,omn_res)
eta = 0.0

AE_HSX = []
AS_HSX = []

for idx, omn in np.ndenumerate(omn_vals):
    file = 'gist_HSX.txt'
    path_to_file = path+file
    gist_class = gist(path_to_file)
    AE_val, AS_val = compute_ae_gist(gist_class,omn,eta*omn)
    AE_HSX         = np.append(AE_HSX,AE_val)
    AS_HSX         = np.append(AS_HSX,AS_val)


plt.loglog(omn_vals,AE_HSX, label='HSX (AE)')
plt.loglog(omn_vals,AS_HSX, label='HSX (AS)')
plt.loglog(omn_vals,AE_HSX-AS_HSX, label='HSX (AE-AS)')
plt.legend()
plt.xlabel('density gradient')
plt.ylabel('available energy')
plt.show()
