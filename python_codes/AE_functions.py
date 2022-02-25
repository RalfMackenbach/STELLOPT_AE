import  numpy           as      np
from    scipy.signal    import  find_peaks
from    numba           import  jit





########################################################################
# This package contains routines to calculate various
# quantities needed for calculating AE
########################################################################




def zerocross1d(x, y, getIndices=True):
  """
    !!! This beautiful little gem is adapted from PyAstronomy !!!
    Find the zero crossing events in a discrete data set.
    Linear interpolation is used to determine the actual
    locations of the zero crossing between two data points
    showing a change in sign. Data point which are zero
    are counted in as zero crossings. Also returns the
    indices left of the zero crossing (even if it hits
    zero exactly).
    *x*
        Array containing the abscissa
    *y*
        Array containing the ordinate
    *getIndices*
        Boolean, if True, also the indicies of the points preceding
        the zero crossing event will be returned. Default is
        True.

    Returns
    ------------------------------------------------------------------
    xvals : array
        The locations of the zero crossing events determined
        by linear interpolation on the data.
    indices : array, optional
        The indices of the points preceding the zero crossing
        events.
  """

  # Check sorting of x-values
  if np.any((x[1:] - x[0:-1]) <= 0.0):
    raise(PE.PyAValError("The x-values must be sorted in ascending order!", \
                         where="zerocross1d", \
                         solution="Sort the data prior to calling zerocross1d."))

  # Indices of points *before* zero-crossing.
  indi = np.where(y[1:]*y[0:-1] < 0.0)[0]

  # Find the zero crossing by linear interpolation
  dx = x[indi+1] - x[indi]
  dy = y[indi+1] - y[indi]
  zc = -y[indi] * (dx/dy) + x[indi]

  # What about the points, which are actually zero?
  # We don't count the first element.
  # Because we are going over a smaller array, the
  # indices are shifted by -1, ensuring we end up
  # left of the crossing in y[zi].
  zi = (np.where(y[1:len(y)] == 0.0)[0])

  # Concatenate indices
  zzindi = np.concatenate((indi, zi))
  # Concatenate zc and locations corresponding to zi
  zz = np.concatenate((zc, x[zi+1]))

  # Sort by x-value
  sind = np.argsort(zz)
  zz, zzindi = zz[sind], zzindi[sind]

  if not getIndices:
    return zz
  else:
    return zz, zzindi



def bounce_wells(l_arr, beta_arr, lambda_var):
    """
        Returns a list of bounce points and indices
        l       = [[-0.1,0.05],[0.31,0.36]]
        l_idx   = [[20,54],[80,94]]
        where each row corresponds to a bounce well
        and the left and right column corresponds to
        the left and right bounce points respectively.
        *l_arr*
            Array of arclength coordinate
        *beta_arr*
            Array of beta values
        *lambda_var*
            pitch-angle
    """
    # Rewrite array in a form suitable for zerocross1d
    zeros_arr   = lambda_var * ( 1 + beta_arr ) - 1

    # Find zero crossings
    l_cross, l_cross_idx = zerocross1d(l_arr, zeros_arr, True)
    # Create empty list in which zero crossings and
    # indices of wells will be stored
    l_wells     = []
    idx_wells   = []

    # Check if zero crossing defines the end or start of a well.
    # If it is a start, it and the next coordinate define a well.
    for (idx, l_idx) in enumerate(l_cross_idx):
        if zeros_arr[l_idx+1]<zeros_arr[l_idx]:
            l_wells.append([l_cross[idx],np.roll(l_cross,-1)[idx]])
            idx_wells.append([l_cross_idx[idx],np.roll(l_cross_idx,-1)[idx]])


    # Return them as numpy arrats
    return l_wells, idx_wells



def integral_I1(l_arr, beta_arr, lambda_var, left_idx, right_idx):
    """
        The integral I1 of all discretely sampled
        points between the bounce points, as
        defined by the indices
        *l_arr*
            Array of arclength coordinate
        *beta_arr*
            Array of beta values
        *lambda_var*
            pitch-angle
        *left_idx*
            the index left of the bounce point,
            such that l[left_idx+1]>=l_bounce
        *bounce_idx*
            The index left of the bounce points so that
            l[right_idx+1]>=l_bounce
    """
    # We define our integration domain for I1
    l_dom    = l_arr[(left_idx+1):(right_idx+1)]
    beta_dom = beta_arr[(left_idx+1):(right_idx+1)]

    # If there is only one point between the bounce points
    # the integral of all interior points vanishes
    if len(l_dom)==1:
        return 0.


    else:
    # Use np roll to shift quantities
        l_dom_shift =    np.roll(l_dom,-1)
        beta_dom_shift = np.roll(beta_dom,-1)

        
        # Delete last index for all
        l_dom           = np.delete(l_dom,-1)
        l_dom_shift     = np.delete(l_dom_shift,-1)
        beta_dom        = np.delete(beta_dom,-1)
        beta_dom_shift  = np.delete(beta_dom_shift,-1)


        # Calculate the integral as obtained by Mathematica
        integral = np.sum((2*(np.sqrt(1 - lambda_var - lambda_var*beta_dom) -  \
        np.sqrt(1 - lambda_var - lambda_var*beta_dom_shift))*(l_dom - \
        l_dom_shift))/(lambda_var*(beta_dom - beta_dom_shift)))
        return integral



def integral_I2(l_arr, beta_arr, f_arr, lambda_var, left_idx, right_idx):
    """
        The integral I2 of all discretely sampled
        points between the bounce points, as
        defined by the indices
        *l_arr*
            Array of arclength coordinate
        *beta_arr*
            Array of beta values
        *f_arr*
            Array of f values
        *lambda_var*
            pitch-angle
        *left_idx*
            the index left of the bounce point,
            such that l[left_idx+1]>=l_bounce
        *bounce_idx*
            The index left of the bounce points so that
            l[right_idx+1]>=l_bounce
    """
    # We define our integration domain for I2
    l_dom    = l_arr[(left_idx+1):(right_idx+1)]
    beta_dom = beta_arr[(left_idx+1):(right_idx+1)]
    f_dom    = f_arr[(left_idx+1):(right_idx+1)]

    # If there is only one point between the bounce points
    # the integral of all interior points vanishes
    if len(l_dom)==1:
        return 0.

    else:
        # Use np roll to shift quantities
        l_dom_shift    = np.roll(l_dom,-1)
        beta_dom_shift = np.roll(beta_dom,-1)
        f_dom_shift    = np.roll(f_dom,-1)

        # Delete last index for all
        l_dom           = np.delete(l_dom,-1)
        l_dom_shift     = np.delete(l_dom_shift,-1)
        beta_dom        = np.delete(beta_dom,-1)
        beta_dom_shift  = np.delete(beta_dom_shift,-1)
        f_dom           = np.delete(f_dom,-1)
        f_dom_shift     = np.delete(f_dom_shift,-1)

        # Calculate the integral as obtained by Mathematica
        integral = np.sum((2*(np.sqrt(1 - lambda_var - lambda_var*beta_dom)*((2 - \
        2*lambda_var + lambda_var*beta_dom - 3*lambda_var*beta_dom_shift)*f_dom \
        + 2*(-1 + lambda_var + lambda_var*beta_dom)*f_dom_shift) + np.sqrt(1 - \
        lambda_var - lambda_var*beta_dom_shift)*(2*(-1 + lambda_var + \
        lambda_var*beta_dom_shift)*f_dom + (2 - 2*lambda_var - \
        3*lambda_var*beta_dom + lambda_var*beta_dom_shift)*f_dom_shift))*(l_dom - \
        l_dom_shift))/(3.*lambda_var**2*(beta_dom - beta_dom_shift)**2))

        return integral



def integral_bound_I1(l_arr, beta_arr, lambda_var, l_b, l_idx, l_0_string):
    """
        The boundary integral of I1
        *l_arr*
            Array of arclength coordinate
        *beta_arr*
            Array of beta values
        *lambda_var*
            pitch-angle
        *l_b*
            the l value of the bounce point
        *l_idx*
            The index left of the bounce points so that
            l[l_idx+1]>=l_b
        *l_0_string*
            The coordinate to which we integrate,
            either "left" or "right"
    """
    # Define relevant quantities
    l_left = l_arr[l_idx]
    l_right= l_arr[l_idx+1]
    beta_left = beta_arr[l_idx]
    beta_right= beta_arr[l_idx+1]
    # Define the integration bound
    if l_0_string == "left":
        integral = (-2*l_b + 2*l_right)/np.sqrt(((beta_left - \
        beta_right)*lambda_var*(l_b - l_right))/(l_left - l_right))
    if l_0_string == "right":
        integral = (2*(l_b - l_left))/np.sqrt(((beta_left - \
        beta_right)*lambda_var*(l_b - l_left))/(l_left - l_right))

    return integral



def integral_bound_I2(l_arr, beta_arr, f_arr, lambda_var, l_b, l_idx, l_0_string):
    """
        The boundary integral of I2
        *l_arr*
            Array of arclength coordinate
        *beta_arr*
            Array of beta values
        *f_arr*
            Array of f values
        *lambda_var*
            pitch-angle
        *l_b*
            the l value of the bounce point
        *l_idx*
            The index left of the bounce points so that
            l[l_idx+1]>=l_b
        *l_0_string*
            The coordinate to which we integrate,
            either "left" or "right"
    """
    # Define relevant quantities
    l_left = l_arr[l_idx]
    l_right= l_arr[l_idx+1]
    beta_left = beta_arr[l_idx]
    beta_right= beta_arr[l_idx+1]
    f_left = f_arr[l_idx]
    f_right= f_arr[l_idx+1]
    # Define the integration bound
    if l_0_string == "left":
        integral = (2.*np.sqrt(((beta_left - beta_right)*lambda_var*(l_b - l_right))/(l_left - l_right))*(2.*f_left*(-l_b + l_right) + f_right*(2.*l_b - 3.*l_left + l_right)))/(3.*(beta_left - beta_right)*lambda_var)
    if l_0_string == "right":
        integral = (2.*(2.*f_right*(-l_b + l_left) + f_left*(2.*l_b + l_left - 3.*l_right))*np.sqrt(((beta_left - beta_right)*lambda_var*(l_b - l_left))/(l_left - l_right)))/(3.*(beta_left - beta_right)*lambda_var)
    return integral




def ground_state(omega_alpha_hat, omega_psi_hat, omega_T_hat,n):
    """
        The ground state equation.
        *omega_alpha_hat*
            normalized bounce averaged drift frequency in the direction of alpha
        *omega_psi_hat*
            normalized bounce averaged drift frequency in the direction of psi
        *omega_T_hat*
            normalized drift wave thingy frequency
        *n*
            the superellipse-parameter, where n>0.
            Also accepts the string "inf" to do the limit n->inf
    """
    if n=="inf":
        numerator   = np.amax(np.asarray([np.abs(omega_T_hat-omega_alpha_hat),np.abs(omega_psi_hat)]))
        denominator = np.amax(np.asarray([np.abs(omega_alpha_hat),np.abs(omega_psi_hat)]))
        result      = numerator/denominator
    elif n>0:
        numerator   = (np.abs(omega_T_hat-omega_alpha_hat)**n + np.abs(omega_psi_hat)**n)**(1/n)
        denominator = (np.abs(omega_alpha_hat)**n + np.abs(omega_psi_hat)**n)**(1/n)
        result      = numerator/denominator
    else:
        raise Exception('n should be > 0, or the string \"inf\". n was {}'.format(n))
    return result



################################################################################
                # SOME FUNCTIONS FOR PERIODIC BOUNDARY CONDITIONS #
################################################################################


def make_periodic_l(l_arr):
    """
        Returns periodic array l_arr_per where,
        l_arr_per   = l_arr.append(l_arr[-1] + Delta_l)
    """
    # We first extend the array to account for the periodic boundary condition
    Delta_l     = np.average( np.delete(np.roll(l_arr,-1),-1) - np.delete(l_arr,-1) )
    l_arr_per   = np.append( l_arr,     l_arr[-1] + Delta_l )

    return l_arr_per






def make_periodic_arr(arr):
    """
        Returns periodic array where,
    """
    # Just latch on first index. This probs doesn't need to be a function
    # but who cares :p
    arr_per=  np.append(arr,  arr[0])

    return arr_per

# wrapper function for periodic boundary conditions on integral I1
# must be fed the arrays after make_periodic runs on them
def I1_periodic(l_arr,beta_arr,lambda_val,l_bounce,l_bounce_idx):
    # Check if the bounce wells cross the periodic boundary

    # If it doesn't, business as usual
    if l_bounce[0]<l_bounce[1]:
        # Calculate integral
        inner_int  =  integral_I1(l_arr, beta_arr, lambda_val,
                        l_bounce_idx[0], l_bounce_idx[1])
        left_bound =  integral_bound_I1(l_arr, beta_arr,
                        lambda_val, l_bounce[0], l_bounce_idx[0], "left")
        right_bound=  integral_bound_I1(l_arr, beta_arr,
                        lambda_val, l_bounce[1], l_bounce_idx[1], "right")
        integral   = inner_int + left_bound + right_bound

    # If it does, we need to be careful about which of the 4 cases we are in
    else:
        # set up handy variable for last l_arr index
        last_l_idx = len(l_arr) - 1
        # Case 1 - the bounce well is fully contained
        if (l_bounce[0]<l_arr[-2]) and (l_bounce[1]>l_arr[1]):
            inner_int_left_per =  integral_I1(l_arr, beta_arr, lambda_val, l_bounce_idx[0], last_l_idx)
            inner_int_right_per=  integral_I1(l_arr, beta_arr, lambda_val, 0, l_bounce_idx[1])
            left_bound =  integral_bound_I1(l_arr, beta_arr,
                            lambda_val, l_bounce[0], l_bounce_idx[0], "left")
            right_bound=  integral_bound_I1(l_arr, beta_arr,
                            lambda_val, l_bounce[1], l_bounce_idx[1], "right")
            integral = inner_int_left_per + inner_int_right_per + left_bound + right_bound
        # Case 2 - the bounce well is left-bounded
        if (l_bounce[0]>l_arr[-2]) and (l_bounce[1]>l_arr[1]):
            inner_int_right_per =  integral_I1(l_arr, beta_arr, lambda_val, 0, l_bounce_idx[1])
            left_bound =  integral_bound_I1(l_arr, beta_arr,
                            lambda_val, l_bounce[0], l_bounce_idx[0], "left")
            right_bound=  integral_bound_I1(l_arr, beta_arr,
                            lambda_val, l_bounce[1], l_bounce_idx[1], "right")
            integral = inner_int_right_per + left_bound + right_bound
        # Case 3 - the bounce well is right-bounded
        if (l_bounce[0]<l_arr[-2]) and (l_bounce[1]<l_arr[1]):
            inner_int_left_per =  integral_I1(l_arr, beta_arr, lambda_val, l_bounce_idx[0], last_l_idx)
            left_bound =  integral_bound_I1(l_arr, beta_arr,
                            lambda_val, l_bounce[0], l_bounce_idx[0], "left")
            right_bound=  integral_bound_I1(l_arr, beta_arr,
                            lambda_val, l_bounce[1], l_bounce_idx[1], "right")
            integral = inner_int_left_per + left_bound + right_bound
        # Case 4 - the bounce well is fully-bounded
        if (l_bounce[0]>l_arr[-2]) and (l_bounce[1]<l_arr[1]):
            left_bound =  integral_bound_I1(l_arr, beta_arr,
                            lambda_val, l_bounce[0], l_bounce_idx[0], "left")
            right_bound=  integral_bound_I1(l_arr, beta_arr,
                            lambda_val, l_bounce[1], l_bounce_idx[1], "right")
            integral = left_bound + right_bound
    return integral




# wrapper function for periodic boundary conditions om integral I2
# must be fed the arrays after make_periodic runs on them
def I2_periodic(l_arr,beta_arr,f_arr,lambda_val,l_bounce,l_bounce_idx):
    # Check if the bounce wells cross the periodic boundary

    # If it doesn't, business as usual
    if l_bounce[0]<l_bounce[1]:
        # Calculate integral
        inner_int  =  integral_I2(l_arr, beta_arr, f_arr, lambda_val,
                        l_bounce_idx[0], l_bounce_idx[1])
        left_bound =  integral_bound_I2(l_arr, beta_arr, f_arr,
                        lambda_val, l_bounce[0], l_bounce_idx[0], "left")
        right_bound=  integral_bound_I2(l_arr, beta_arr, f_arr,
                        lambda_val, l_bounce[1], l_bounce_idx[1], "right")
        integral   = inner_int + left_bound + right_bound

    # If it does, we need to be careful about which of the 4 cases we are in
    else:
        # set up handy variable of last l_arr index
        last_l_idx = len(l_arr) - 1
        # Case 1 - the bounce well is fully contained
        if (l_bounce[0]<l_arr[-2]) and (l_bounce[1]>l_arr[1]):
            inner_int_left_per =  integral_I2(l_arr, beta_arr, f_arr, lambda_val, l_bounce_idx[0], last_l_idx)
            inner_int_right_per=  integral_I2(l_arr, beta_arr, f_arr, lambda_val, 0, l_bounce_idx[1])
            left_bound =  integral_bound_I2(l_arr, beta_arr, f_arr,
                            lambda_val, l_bounce[0], l_bounce_idx[0], "left")
            right_bound=  integral_bound_I2(l_arr, beta_arr, f_arr,
                            lambda_val, l_bounce[1], l_bounce_idx[1], "right")
            integral = inner_int_left_per + inner_int_right_per + left_bound + right_bound
        # Case 2 - the bounce well is left-bounded
        if (l_bounce[0]>l_arr[-2]) and (l_bounce[1]>l_arr[1]):
            inner_int_right_per =  integral_I2(l_arr, beta_arr, f_arr, lambda_val, 0, l_bounce_idx[1])
            left_bound =  integral_bound_I2(l_arr, beta_arr, f_arr,
                            lambda_val, l_bounce[0], l_bounce_idx[0], "left")
            right_bound=  integral_bound_I2(l_arr, beta_arr, f_arr,
                            lambda_val, l_bounce[1], l_bounce_idx[1], "right")
            integral = inner_int_right_per + left_bound + right_bound
        # Case 3 - the bounce well is right-bounded
        if (l_bounce[0]<l_arr[-2]) and (l_bounce[1]<l_arr[1]):
            inner_int_left_per =  integral_I2(l_arr, beta_arr, f_arr, lambda_val, l_bounce_idx[0], last_l_idx)
            left_bound =  integral_bound_I2(l_arr, beta_arr, f_arr,
                            lambda_val, l_bounce[0], l_bounce_idx[0], "left")
            right_bound=  integral_bound_I2(l_arr, beta_arr, f_arr,
                            lambda_val, l_bounce[1], l_bounce_idx[1], "right")
            integral = inner_int_left_per + left_bound + right_bound
        # Case 4 - the bounce well is fully-bounded
        if (l_bounce[0]>l_arr[-2]) and (l_bounce[1]<l_arr[1]):
            left_bound =  integral_bound_I2(l_arr, beta_arr, f_arr,
                            lambda_val, l_bounce[0], l_bounce_idx[0], "left")
            right_bound=  integral_bound_I2(l_arr, beta_arr, f_arr,
                            lambda_val, l_bounce[1], l_bounce_idx[1], "right")
            integral = left_bound + right_bound
    return integral







################################################################################
                # SOME FUNCTIONS NEEDED IN THE CALCULATION #
###############################################################################

def lambda_filtered(lambda_arr,beta_arr,delta_lambda):
    """
        Returns filtered array, where all values in
        lambda arr which lie within delta_lambda*range(beta)
        are deleted.

        lambda_arr      - array with lambda values
        beta_arr        - array with beta values
        delta_lambda    - the singularity padding
    """
    #make new container
    lambda_arr_new  = lambda_arr
    # find idx of local maxima of beta arrays
    beta_max_idx    = find_peaks(beta_arr)[0]
    # Find corresponding beta vals and lambda vals
    beta_local_max  = np.asarray([beta_arr[i] for i in beta_max_idx])
    lamdba_inf      = 1.0 / ( 1.0 + beta_local_max )
    # construct range(lambda)
    max_beta    = np.amax(beta_arr)
    min_beta    = np.amin(beta_arr)
    lambda_max  = 1.0/(1.0+min_beta)
    lambda_min  = 1.0/(1.0+max_beta)
    lambda_range= lambda_max-lambda_min
    # loop over lambda inf and delete lambdas within singularity padding
    for lambda_inf_val in lamdba_inf:
        lower_bound         =   lambda_inf_val - delta_lambda*lambda_range
        upper_bound         =   lambda_inf_val + delta_lambda*lambda_range
        lambda_delete_idx   =   np.where(np.logical_and(lambda_arr_new>=lower_bound, lambda_arr_new<=upper_bound))
        lambda_arr_new      =   np.delete(lambda_arr_new, lambda_delete_idx)

    return lambda_arr_new



def normalized_frequencies(l_arr,beta_arr,lambda_arr,dbetadpsi,dbetadalpha,Delta_psi,Delta_alpha):
    """
        Returns the normalized bounce frequency and bounce time
        omega_hat_alpha, omega_hat_psi, J_hat.
        These are lists of lists where the first list index
        corresponds to the lambda coordinate, and the second
        index denotes the bounce well.
        Must be fed arrays after make_periodic runs on them.

        l_arr           - array with arclength values
        beta_arr        - array with beta values
        lambda_arr      - array with lambda values
        dbetadpsi       - array with dbetadpsi values
        dbetadalpha     - array with dbetadalpha values
        Delta_psi       - box size in psi direction
        Delta_alpha     - box size in alpha direction
    """
    # Now let's calculate the normalized bounce frequences for
    # all the different lambda values.
    # We store all these in a list of lists, where the primary
    # list idx corresponds to the lambda_arr_idx, and the secondary idx
    # corresponds to the bounce well.
    omega_alpha = []
    omega_psi   = []
    J_hat       = []
    L_total     = l_arr[-1]-l_arr[0]
    # Loop over all lambda vals
    for lambda_val in lambda_arr:
        # Create containers in which to store each well and
        # calculate the bounce points
        omega_alpha_container = []
        omega_psi_container   = []
        J_hat_container       = []
        l_bounce_wells, l_bounce_indices   = bounce_wells(l_arr, beta_arr, lambda_val)
        # Loop over each bounce well
        for (l_idx, l_bounce) in enumerate(l_bounce_wells):
            l_bounce_idx = l_bounce_indices[l_idx]
            omega_alpha_numerator = I2_periodic(l_arr,beta_arr,dbetadpsi,lambda_val,l_bounce,l_bounce_idx)
            omega_psi_numerator   = I2_periodic(l_arr,beta_arr,dbetadalpha,lambda_val,l_bounce,l_bounce_idx)
            denominator           = I1_periodic(l_arr,beta_arr,lambda_val,l_bounce,l_bounce_idx)
            # Calculate the bounce averaged frequencies and add to containers
            omega_alpha_val           = lambda_val * Delta_psi * omega_alpha_numerator/denominator
            omega_alpha_container.append(omega_alpha_val)
            omega_psi_val             = - lambda_val * Delta_alpha * omega_psi_numerator/denominator
            omega_psi_container.append(omega_psi_val)
            # Calculate J_hat and add to container
            J_hat_val                 = denominator/L_total
            J_hat_container.append(J_hat_val)
        # Append containers to lists
        omega_alpha.append(omega_alpha_container)
        omega_psi.append(omega_psi_container)
        J_hat.append(J_hat_container)

    return omega_alpha, omega_psi, J_hat






def construct_integrand(z_arr,lambda_arr,omega_alpha,omega_psi,J_hat,grad_T,grad_n,Delta_psi,Delta_alpha,n_geom):
    # Create multidimensional grid
    gridpoints   = np.meshgrid(z_arr,lambda_arr)
    z_grid       = gridpoints[0]
    lambda_grid  = gridpoints[1]
    # Create array holding the integrand
    integrand    = np.empty_like(z_grid)

    # Iterate through all indices
    for index, x in np.ndenumerate(z_grid):
        lambda_idx = index[0]
        z_val      = z_grid[index]
        lambda_val = lambda_grid[index]
        integrand_wells  = []
        # Calculate integrand per well
        for idx, omega_alpha_well in enumerate(omega_alpha[lambda_idx]):
            omega_alpha_well = (omega_alpha[lambda_idx])[idx]
            omega_psi_well   = (omega_psi[lambda_idx])[idx]
            J_hat_well       = (J_hat[lambda_idx])[idx]
            omega_T_hat      = Delta_psi*(1.0 * grad_T / z_val + grad_n * (1.0 - 3.0/2.0 * 1.0/z_val) )
            ground_state_val = ground_state(omega_alpha_well, omega_psi_well, omega_T_hat,n_geom)
            integrand_val    = np.exp(-z_val)*z_val**(2.5) * ( omega_alpha_well**2.0 *
                              ( omega_T_hat/omega_alpha_well - 1 + ground_state_val ) +
                              omega_psi_well**(2.0) * (-1 + ground_state_val) ) * J_hat_well
            integrand_wells.append(integrand_val)
        # Sum all wells
        integrand[index] = np.sum(integrand_wells)

    return z_grid, lambda_grid, integrand
