# Barion Number Violation Running (bnv-running)

This project numerically solves the coupled differential equations for Wilson coefficients and gauge/Yukawa couplings in the context of effective field theory (EFT) scenarios. It uses Python and `scipy`'s ODE solver to evolve parameters across energy scales.

---

## 📁 Project Structure
```
├── results/ # (Optional) Folder to store output results or figures  
├── wce/  
│ ├── __init__.py # Package initializer  
│ └── wce.py # Main module for evolution of couplings and Wilson coefficients  
├── example.py # Example script demonstrating usage  
├── n-n_bar-oscillation.ipynb # Use case: Neutron–antineutron oscillation example  
└── requirements.txt # Python dependencies
```
---

## ⚙️ Features

- Solves renormalization group equations (RGEs) for:  
  - Yukawa couplings: `yu`, `yd`  
  - Gauge couplings: `g1`, `g2`, `g3`  
  - 9 Wilson coefficients  
- Piecewise evolution across multiple energy regions with support for model structure changes  
- Easy integration into research pipelines using Jupyter Notebooks or scripts  
- Modular and extensible code structure

---

## 🛠️ Requirements

Install all dependencies with:

```bash
pip install -r requirements.txt
```

## 🚀 How to Use

### 🔄 ```evolve(...)``` Function – Inputs and Outputs

The core of this project lies in the evolve function from wce/wce.py. This function numerically solves the coupled differential equations for Wilson coefficients and couplings across multiple energy regions.

#### Inputs to ```evolve(...)```:
```
t_solutions, y_solutions, all_results = evolve(E_transitions, NJKL_val, X_val, initial_couplings, initial_wc, E_END)
```

- E_transitions:  
A list of energies (in eV) where new physics is introduced — at each energy boundary, the operator or model structure can change.

- NJKL_val:  
A list of [N, J, K, L] values at each energy region. These integers define the active structure of operators or couplings.

Set an index to -1 to remove its contribution in that region.

- X_val:  
List indicating whether certain special operators (like X-type) are active (e.g., 0 or 1) in each region.

- initial_couplings:  
A nested list of coupling values (yu, yd, g1, g2, g3) at each energy transition. (Read below for detail explanation!)

Use None to inherit from the previous energy region automatically.

- initial_wc:  
A list of 9 floating-point values specifying the initial Wilson coefficients at the starting energy.

- E_END:  
Final energy (in eV) up to which the RG evolution should proceed.


#### Outputs from ```evolve(...)```:
- t_solutions:  
List of 1D NumPy arrays, each representing the log(E) time grid for an energy region.

- y_solutions:  
List of dictionaries with evolved quantities in each region. Each dictionary contains:

- 'wc':  
Wilson coefficients

- 'yu', 'yd':  
Yukawa couplings

- 'g1', 'g2', 'g3':  
Gauge couplings

- all_results:  
Full internal output of the solver, including all numerical states — useful for advanced analysis.


### 🧩 Understanding ```initial_couplings```

The ```initial_couplings``` input defines the values of Yukawa and gauge couplings at each energy region.

Each item in ```initial_couplings``` corresponds to the starting values at a specific energy transition point, and contains values for: 

```
[
  [y_u],
  [y_d],
  [g1],
  [g2],
  [g3]
]
```

Each of these is a list of numbers — the length of the list represents how many couplings of that type are active in that region.
Each value inside the list is either a number (explicit initial value) or None (to inherit from previous evolution).


#### 🔄 How it works:

Let’s say for a given energy region we have:

```
initial_couplings_1 = [
    [None, 1],  # y_u: keep the evolved value for y_u0, set y_u1 = 1
    [None, 1],  # y_d: same logic
    [None, 1],  # g1: keep g_10, set g_11 = 1
    [None],     # g2: keep g_20
    [None, 1],  # g3: keep g_30, set g_31 = 1
]
```
Here’s how to interpret this:

- At energy transition E[1], y_u has 2 components: y_u0, y_u1

      - y_u0 is inherited from the evolved result at E[0]

      - y_u1 is newly introduced and initialized to 1

- g1 now has two active components: g_10, g_11

      - g_10: inherited

      - g_11: explicitly set to 1

- For any position where you write None, it uses the result from the previous energy range evolution.

#### 🔢 General Rules

- Each inner list (e.g. [None, 1]) corresponds to that coupling's components in the current energy region.

- The position in the list (0th, 1st, 2nd, …) indicates the index of the coupling:  
e.g., g_3l → g_30, g_31, g_32, etc.

- None = "inherit evolved value from previous region"

- Number = "explicitly set this value at the start of this region"
  

### 📌 Check the Examples

Refer to ```example.py``` for a complete working demonstration of how to use the package.

It covers:

- Setting energy transition points at which model structure changes (new operators turn on/off)

- Specifying initial values and matching conditions for Yukawa/gauge couplings across energy scales

- Evolving Wilson coefficients and plotting their energy dependence
