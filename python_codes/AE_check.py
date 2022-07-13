import  f90nml
import  AE_routines     as  ae
import  re
import  numpy as np
import  matplotlib.pyplot   as plt

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
    lam_res = 1000
    Delta_theta = 1e-10
    del_sing = 0.0
    ae_val = ae.ae_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr,lam_res,Delta_theta,del_sing,L_tot,omnigenous)
    return ae_val




################################ DO EXAMPLE CASE ###############################
# If one expands AE in a weakly driven regime, one finds a square dependence on
# the gradients. In a strongly driven regime, one expects a linear dependence.
# Check if that is what you find! This picture changes for devices which have
# omega_psi = 0

path = 'gist_files/'
omn_res=100
omn_hlf=int(omn_res/2)
omn_vals = np.logspace(-6,6,omn_res)
eta = 1.0
AE_W7X = []
AE_D3D = []
AE_HSX = []
AE_NCSX= []

for idx, omn in np.ndenumerate(omn_vals):
    file = 'gist_W7XSC.txt'
    path_to_file = path+file
    gist_class = gist(path_to_file)
    AE_val = compute_ae_gist(gist_class,omn,eta*omn,omnigenous=False)
    AE_W7X = np.append(AE_W7X,AE_val)

for idx, omn in np.ndenumerate(omn_vals):
    file = 'gist_D3D.txt'
    path_to_file = path+file
    gist_class = gist(path_to_file)
    AE_val = compute_ae_gist(gist_class,omn,eta*omn,omnigenous=True)
    AE_D3D = np.append(AE_D3D,AE_val)

for idx, omn in np.ndenumerate(omn_vals):
    file = 'gist_HSX.txt'
    path_to_file = path+file
    gist_class = gist(path_to_file)
    AE_val = compute_ae_gist(gist_class,omn,eta*omn,omnigenous=False)
    AE_HSX = np.append(AE_HSX,AE_val)

for idx, omn in np.ndenumerate(omn_vals):
    file = 'gist_NCSX.txt'
    path_to_file = path+file
    gist_class = gist(path_to_file)
    AE_val = compute_ae_gist(gist_class,omn,eta*omn,omnigenous=False)
    AE_NCSX = np.append(AE_NCSX,AE_val)


plt.loglog(omn_vals,AE_W7X, label='W7X')
plt.loglog(omn_vals,AE_D3D, label='D3D')
plt.loglog(omn_vals,AE_HSX, label='HSX')
plt.loglog(omn_vals,AE_NCSX,label='NCSX')
plt.loglog(omn_vals[(omn_hlf)::],omn_vals[(omn_hlf)::], label='1-law',linestyle='dashdot',color='black')
plt.loglog(omn_vals[0:(omn_hlf)],omn_vals[0:(omn_hlf)]**(2), label='2-law',linestyle='dashed',color='black')
plt.loglog(omn_vals[0:int((omn_hlf/2))],omn_vals[0:int((omn_hlf/2))]**(7/2)*10, label='7/2-law',linestyle='dotted',color='black')
plt.legend()
plt.xlabel('omn')
plt.ylabel('available energy')
plt.show()
