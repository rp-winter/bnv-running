import numpy as np
from scipy.integrate import solve_ivp 
import math

# DE for the couplings
def couplings_de(yu, yd, g1, g2, g3, NJKL, X):
    Nf = 6 # flavor number
    N = NJKL[0]
    J = NJKL[1]
    K = NJKL[2]
    L = NJKL[3]
    
    # beta functions for the couplings
    sum_g1_sq = np.sum(g1**2)
    sum_g2_sq = np.sum(g2**2)
    sum_g3_sq = np.sum(g3**2)

    sum_yu_sq = np.sum(yu[1:]**2)
    sum_yd_sq = np.sum(yd[1:]**2)

    Yu = np.array([
        [1.3e-5, 0, 0],
        [0, 7.3e-3, 0],
        [0, 0, 1.0]
    ])

    Yd = np.array([
        [2.7e-5, 0, 0],
        [0, 5.5e-4, 0],
        [0, 0, 0.024]
    ])

    T = np.trace(3*np.matmul(Yd.T, Yd) + 3*np.matmul(Yu.T, Yu))

    beta1 = (-(2/3)*Nf - (N+1)/10)*g1**2
    beta2 = g2**2 * (16/3 - (2/3)*Nf - (N+1)/6) + 2*sum_g2_sq
    beta3 = g3**2 * (8 - 2*Nf/3) + 3*sum_g3_sq

    beta_yu = (3/2)*(Yu[0, 0]**2 - Yd[0, 0]**2) + T + ((9/2)*sum_yu_sq + (3/2)*sum_yd_sq) - ( (17/20)*sum_g1_sq + (9/4)*sum_g2_sq + 8*sum_g3_sq )
    beta_yd = (3/2)*(Yd[0, 0]**2 - Yu[0, 0]**2) + T + ((9/2)*sum_yd_sq + (3/2)*sum_yu_sq) - ( (1/4)*sum_g1_sq + (9/4)*sum_g2_sq + 8*sum_g3_sq )

    # return the derivatives
    # yu0, yu1, ..., yuN, yd0, yd1, ..., ydN, g10, g11, ..., g1J, g20, g21, ..., g2K, g30, g31, ..., g3L
    dcouplings_dt = []
    sixteen_pi_sq = 16*math.pi**2
    for i in range(N+1):
        dcouplings_dt.append(beta_yu*yu[i]/sixteen_pi_sq)
    for i in range(N+1):
        dcouplings_dt.append(beta_yd*yd[i]/sixteen_pi_sq)
    dcouplings_dt.extend(-beta1*g1/sixteen_pi_sq)
    dcouplings_dt.extend(-beta2*g2/sixteen_pi_sq)
    dcouplings_dt.extend(-beta3*g3/sixteen_pi_sq)
    return dcouplings_dt

# DE for the Wilson coefficients
def wc_de(yu, yd, g1, g2, g3, wc, NJKL):
    N = NJKL[0]
    J = NJKL[1]
    K = NJKL[2]
    L = NJKL[3]
    
    sum_g1_sq = np.sum(g1**2)
    sum_g2_sq = np.sum(g2**2)
    sum_g3_sq = np.sum(g3**2)

    sum_yu_sq = np.sum(yu**2)
    sum_yd_sq = np.sum(yd**2)

    sum_yu_yd = np.sum(yu*yd)

    c1 = 3*(sum_yu_sq + 2*sum_yd_sq)/2

    # Diagonal elements {The (4 pi)^2 are factored out}
    gamma11 = -sum_g1_sq
    gamma22 = (-sum_g1_sq + 9*sum_g2_sq + 8*sum_g3_sq + c1 )
    gamma33 = ( sum_g1_sq / 2 + 9*sum_g2_sq/2 + 8*sum_g3_sq + c1 )
    gamma44 = (-4*sum_g1_sq + c1)
    gamma55 = (2*sum_g1_sq + c1)
    gamma66 = ( 2*sum_g1_sq + 8*sum_g3_sq + c1 )
    gamma77 = (-sum_g1_sq + 9*sum_g2_sq + 8*sum_g3_sq + c1)
    gamma88 = ( sum_g1_sq / 2 + 9*sum_g2_sq/2 + 8*sum_g3_sq + c1)
    gamma99 = (2*sum_g1_sq + c1)

    # Off-diagonal elements
    gamma17 = gamma32 = gamma87 = 6 * sum_yu_yd
    gamma23 = gamma58 = gamma78 = gamma93 = 12 * sum_yu_yd
    gamma48 = gamma63 = 18 * sum_yu_yd
    gamma39 = gamma85 = 3 * sum_yu_yd

    Gamma = [
        [gamma11, 0, 0, 0, 0, 0, gamma17, 0, 0],
        [0, gamma22, gamma23, 0, 0, 0, 0, 0, 0],
        [0, gamma32, gamma33, 0, 0, 0, 0, 0, gamma39],
        [0, 0, 0, gamma44, 0, 0, 0, gamma48, 0],
        [0, 0, 0, 0, gamma55, 0, 0, gamma58, 0],
        [0, 0, gamma63, 0, 0, gamma66, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, gamma77, gamma78, 0],
        [0, 0, 0, 0, gamma85, 0, gamma87, gamma88, 0],
        [0, 0, gamma93, 0, 0, 0, 0, 0, gamma99],
    ]

    # Return the derivatives
    return np.dot(np.transpose(Gamma), wc)/(16*math.pi**2)
    
def seperate_states(state, N, J, K, L, O):
    """
    Seperates the state vector into yu, yd, g1, g2, g3 and wilson coefficients.
    :param state: state vector
    :param N: N
    :param J: J
    :param K: K
    :param L: L
    :param O: number of wilson coefficients
    """
    # check validity of state vector length
    if len(state) != 2*(N+1) + (J+1) + (K+1) + (L+1) + O:
        raise ValueError("Invalid state vector length. Expected length: {}, got: {}".format(2*(N+1) + (J+1) + (K+1) + (L+1) + O, len(state)))
    # Seperate the states into yu, yd, g1, g2, g3 and wilson coefficients
    # Extract variables from state vector
    yu = state[:(N+1)]       # First N elements
    yd = state[(N+1):2*(N+1)]    # Next N elements
    g1 = state[2*(N+1):2*(N+1)+(J+1)]  # Next J elements
    g2 = state[2*(N+1)+(J+1):2*(N+1)+(J+1)+(K+1)]  # Next K elements
    g3 = state[2*(N+1)+(J+1)+(K+1):2*(N+1)+(J+1)+(K+1)+(L+1)]  # Next L elements
    wc = state[-O:]  # Last O elements for wilson coefficients

    return yu, yd, g1, g2, g3, wc

# Combined ODE function for solve_ivp
def combined_odes(t, state, NJKL, X):
    """"
    t : float
        log(Energy)
        Not used in the function but required by solve_ivp
    state : array_like
        State vector containing couplings and wilson coefficients
    NJKL : list
        List containing N, J, K, L values
    X: Value of X at the current energy
    """
    N, J, K, L = NJKL

    yu, yd, g1, g2, g3, wc = seperate_states(state, N, J, K, L, 9)

    # Compute derivatives
    d_couplings = couplings_de(yu, yd, g1, g2, g3, NJKL, X) # Returns array of length (2N + J + K + L)
    d_wc = wc_de(yu, yd, g1, g2, g3, wc, NJKL) # Matrix-vector product for 9 elements

    # Return concatenated derivatives
    return np.concatenate([d_couplings, d_wc])

def evolve(E_transitions, NJKL_val, X_val, initial_couplings, initial_wc, E_END, grid_points=500):
    """
    Solves the coupled ODEs for couplings & wilson coefficients across energy regions.
    
    Parameters:
    - E_transitions: List of transition energy points.
    - NJKL_val: List of (N, J, K, L) values for each energy region.
    - initial_couplings: Initial values for couplings.
    - initial_wc: Initial values for wilson coefficients.
    - E_END: Final energy scale.
    - grid_points: Either an integer (same grid size for all regions) or a list (custom grid size per region).
    """

    # Convert Energy to log scale
    t_transitions = np.log(E_transitions)
    t_final = np.log(E_END)

    # Store solutions
    t_solutions = []
    y_solutions = []
    sol_info = []

    # Ensure grid_points is a list for per-region customization
    if isinstance(grid_points, int):
        grid_points = [grid_points] * (len(E_transitions))  # Use same value for all regions
    elif len(grid_points) != len(E_transitions):
        raise ValueError("grid_points list length must match the number of energy regions.")

    # Initial conditions (couplings + wilson coefficients)
    y_start = np.concatenate(initial_couplings[0])  # Flatten initial couplings
    y_start = np.concatenate((y_start, initial_wc))  # Add wilson coefficients

    # Start at the first energy point
    t_start = t_transitions[0]

    O = len(initial_wc)  # Number of wilson coefficients

    # Solve in each region
    for i in range(len(E_transitions) - 1):
        t_end = t_transitions[i + 1]
        NJKL = NJKL_val[i]
        N, J, K, L = NJKL
        X = X_val[i]

        # Fine grid for evaluation (user-defined, either uniform or per-region)
        t_fine_eval = np.linspace(t_start, t_end, grid_points[i])

        # Solve the system in this region
        sol = solve_ivp(combined_odes, [t_start, t_end], y_start, args=(NJKL, X), t_eval=t_fine_eval)
        
        # Store results
        t_solutions.append(sol.t)
        y_solutions.append(sol.y)
        sol_info.append({
            'status': sol.status,
            'success': sol.success,
            'message': sol.message,
        })

        # Update initial conditions for next region
        if (i + 1) < len(E_transitions):  
            y_last = sol.y[:, -1]  # Get final values from this region
            new_initial_conditions = []

            # Update couplings if provided
            yu_last, yd_last, g1_last, g2_last, g3_last, wc = seperate_states(y_last, N, J, K, L, O)

            new_NJKL = NJKL_val[i + 1]
            new_NJKL_extended = [new_NJKL[0], new_NJKL[0], new_NJKL[1], new_NJKL[2], new_NJKL[3]]
            coupling_types = [yu_last, yd_last, g1_last, g2_last, g3_last]

            for idx, last_vals in enumerate(coupling_types):
                for j in range(new_NJKL_extended[idx]+1):
                    new_val = initial_couplings[i + 1][idx][j]
                    new_initial_conditions.append(last_vals[j] if new_val is None else new_val)
           
            # keep wilson coefficients evolving
            y_start = np.concatenate((new_initial_conditions, wc))  # Combine new couplings with last wilson coefficients
        t_start = t_end  # Move to next region

    # Solve the final region from last transition point to E_END
    NJKL = NJKL_val[-1]
    X = X_val[-1]
    t_fine_eval = np.linspace(t_start, t_final, grid_points[-1])  # Use last grid setting
    sol = solve_ivp(combined_odes, [t_start, t_final], y_start, args=(NJKL,X), t_eval=t_fine_eval)

    t_solutions.append(sol.t)
    y_solutions.append(sol.y)
    sol_info.append({
            'status': sol.status,
            'success': sol.success,
            'message': sol.message,
        })

    # make each element of y_solution a dictionary with the corresponding counplings and wilson coefficeints seperated
    for i in range(len(y_solutions)):
        yu, yd, g1, g2, g3, wc = seperate_states(y_solutions[i], NJKL_val[i][0], NJKL_val[i][1], NJKL_val[i][2], NJKL_val[i][3], O)
        y_solutions[i] = {
            'yu': yu,
            'yd': yd,
            'g1': g1,
            'g2': g2,
            'g3': g3,
            'wc': wc
        }

    return t_solutions, y_solutions, sol_info
