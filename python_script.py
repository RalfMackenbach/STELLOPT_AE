import  AE_functions        as      AE
import  numpy               as      np
from    scipy.special       import  gamma
import  G_functions         as      pyG
from    scipy.integrate     import  cumtrapz
from    scipy.interpolate   import  interp1d


def calculate_AE(B_arr, dBdx_arr, dBdy_arr, l_arr, Delta_x, Delta_y, omn, omt,
                 lambda_res, delta_lambda, z_res, z_min, z_max, n_geom=2):
    """
    Calculates available energy of a flux-tube, using GENE classes

    Inputs in order:
    B_arr       - array of B as a function of arclength
    dBdx_arr    - array of dBdx as a function of arclength
    dBdy_arr    - array of dBdy as a function of arclength
    l_arr       - array of arclengths
    Delta_x     - width of fluxtube in x
    Delta_y     - width of fluxtube in y
    omn         - density gradient, same as GENE
    omt         - temperature gradient, same as GENE
    lambda_res  - resolution of lambda grid
    delta_lambda- distance between singularties
    z_res       - resolution of z_grid
    z_min       - minimal z value of z_grid
    z_max       - maximal z value of z_grid
    n_geom      - geometry parameter of flux tube (superellipse),
                  accepts \"inf\" for square flux tube. Should be 2,
                  as this is the only valid solution.

    returns the dimensionless AE
    """

    # Define gradients
    grad_T          = -omt
    grad_n          = -omn

    # Rescale arrays
    L_tot           = l_arr[-1]-l_arr[0]
    B_ave           = np.trapz(B_arr, x=l_arr)/L_tot
    beta            = B_arr/B_ave
    beta_x          = dBdx_arr/B_ave
    beta_y          = dBdy_arr/B_ave

    # Make arrays periodic
    beta_per        = np.append(beta,  beta[0])
    dbetadx_per     = np.append(beta_x,beta_x[0])
    dbetady_per     = np.append(beta_y,beta_y[0])
    Delta_l         = np.average( np.delete(np.roll(l_arr,-1),-1) - np.delete(l_arr,-1) )
    l_hat_per       = np.append( l_arr,     l_arr[-1] + Delta_l )/L_tot

    # Make lambda array
    max_beta        = np.amax(beta_per)
    min_beta        = np.amin(beta_per)
    lambda_max      = 1.0/(1.0+min_beta)
    lambda_min      = 1.0/(1.0+max_beta)
    lambda_range    = lambda_max-lambda_min
    lambda_arr      = np.linspace(lambda_min,lambda_max,num=int(lambda_res),
                                  endpoint=False)

    # Filter singularities out of lambda array
    lambda_arr_filt = AE.lambda_filtered(lambda_arr,beta_per,delta_lambda)

    # Now let's calculate the normalized bounce frequences for
    # all the different lambda values.
    # We store all these in a list of lists, where the primary
    # list idx corresponds to the lambda_arr idx, and the secondary idx
    # corresponds to the bounce well.
    omega_alpha, omega_psi, J_hat = AE.normalized_frequencies(l_hat_per,beta_per,
                                    lambda_arr_filt,dbetadx_per,
                                    dbetady_per,Delta_x,Delta_y)

    # Let's construct our integrand
    # We create a z_array for our dimensionless energy parameter
    z_arr = np.linspace(z_min,z_max,num=int(z_res),endpoint=True)
    z_grid, lamdba_grid, integrand = AE.construct_integrand(z_arr,
                                     lambda_arr_filt,omega_alpha,omega_psi,
                                     J_hat,grad_T,grad_n,Delta_x,Delta_y,
                                     n_geom)

    # And finally do trapz twice to find AE
    integral_over_z     = np.trapz(np.asarray(integrand),np.asarray(z_grid),axis=-1)
    integral_full       = np.trapz(integral_over_z,np.asarray(lambda_arr_filt))

    dimless_AE          = integral_full

    return dimless_AE
