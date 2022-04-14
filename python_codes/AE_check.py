import  f90nml
import  AE_routines     as  ae
import  re
import  numpy as np

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
    dBdx        = gist_class.functions[:,6]
    dBdy        = gist_class.functions[:,5]
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
    z_res = 1000
    lam_res = 1000
    Delta_theta = 1e-10
    return ae.ae_total(q0,dlnTdx,dlnndx,Delta_x,Delta_y,b_arr,dbdx_arr,dbdy_arr,sqrt_g,theta_arr,lam_res,z_res,z_min,z_max,Delta_theta,L_tot)



# import all gist classes
path = '/Users/ralfmackenbach/Documents/GitHub/STELLOPT_AE/python_codes/gist_files/'
files = sorted(['gist_D3D.txt','gist_NCSX.txt','gist_W7XHM.txt','gist_HSX.txt','gist_STELLOPT.txt','gist_W7XSC.txt'])

for file in files:
    path_to_file = path+file
    gist_class = gist(path_to_file)
    print(compute_ae_gist(gist_class,4.0,0.0))
