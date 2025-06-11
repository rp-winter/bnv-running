# Barion Number Violation Running (bnv-running)

This project numerically solves the coupled differential equations for Wilson coefficients and gauge/Yukawa couplings in the context of effective field theory (EFT) scenarios. It uses Python and `scipy`'s ODE solver to evolve parameters across energy scales.

---

## 📁 Project Structure
```
├── bnv_running/ 
│   ├── bnv_running.py
│   ├── differential_equations.py 
├── n-n_bar-oscillation
│   ├── results/
│   ├── plot.py
│   ├── wc_evol_and_lambda.ipynb
└── requirements.txt # Python dependencies
```
---

## ⚙️ Features

- Solves renormalization group equations (RGEs) for:  
  - Yukawa couplings: `yu`, `yd`  
  - Gauge couplings: `g1`, `g2`, `g3`  
  - Wilson coefficients  
- Piecewise evolution across multiple energy regions with support for model structure changes  

---

## 🚀 Getting Started
**Clone the repository**:
   ```bash
   git clone https://github.com/rp-winter/bnv-running.git
   cd bnv-running
   ```

### 🛠️ Requirements

Install all dependencies with:

```bash
pip install -r requirements.txt
```

## 🚀 How to Use

A full explanation of all modules, functions, and physics context is available in the [Wiki](https://github.com/rp-winter/bnv-running/wiki)


Refer to ```n-n_bar-oscillation/wc_evol_and_lambda.ipynb``` for a working example.
