import numpy as np
from scipy.integrate import solve_ivp 
from differential_equations import couplings_de, wc_de

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
def combined_odes(t, state, NJKL, gamma_matrix_fn, wc_de_kwargs):

    N, J, K, L = NJKL

    yu, yd, g1, g2, g3, wc = seperate_states(state, N, J, K, L, 9)

    # Compute derivatives
    d_couplings = couplings_de(yu, yd, g1, g2, g3, NJKL) # Returns array of length (2N + J + K + L)
    d_wc = wc_de(t, wc, yu, yd, g1, g2, g3, NJKL, gamma_matrix_fn, wc_de_kwargs)  # Returns array of length O

    # Return concatenated derivatives
    return np.concatenate([d_couplings, d_wc])

def evolve(E_transitions, NJKL_val, initial_couplings, initial_wc, E_END, gamma_matrix_fn, grid_points=500, wc_de_kwargs={}):
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

        # Fine grid for evaluation (user-defined, either uniform or per-region)
        t_fine_eval = np.linspace(t_start, t_end, grid_points[i])

        # Solve the system in this region
        def ode_system(t, y):
            return combined_odes(t, y, NJKL, gamma_matrix_fn, wc_de_kwargs)
        sol = solve_ivp(ode_system, [t_start, t_end], y_start, t_eval=t_fine_eval)
        
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
    t_fine_eval = np.linspace(t_start, t_final, grid_points[-1])  # Use last grid setting

    def ode_system(t, y):
            return combined_odes(t, y, NJKL, gamma_matrix_fn, wc_de_kwargs)
    sol = solve_ivp(ode_system, [t_start, t_final], y_start, t_eval=t_fine_eval)

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
