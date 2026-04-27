# DyAux - Dynamic Design and Optimization of Auxotrophic Communities
This algorithm helps build auxotrophic communities and further optimize them for bioproduction

---

## Key Scripts  
- **Masterscripts for Four Environments (effect of medium composition and oxygen availability):**
  - `AuxMasterHetero.m`  
  - `AuxMasterHomo.m`  
  - `FSEOF_dfba_masterscript.m`  
  
   Initial inputs like models and parameters (environment and kinetic) can be changed by manipulating these files.

---

## Models Used  
All models were obtained from the **BiGG Models** and **BioModels** databases. Some annotations were adjusted to ensure consistent formatting across the dataset.

---

## Prerequisites  
- **MATLAB**
- **COBRA Toolbox**  
- **IBM ILOG CPLEX**

All simulations were conducted using **MATLAB R2018a**, **COBRA Toolbox v3.0**, and **IBM ILOG CPLEX 12.8**.

Please initialise the COBRA Toolbox and the solver before running DyAux
We recommend the ibm_cplex solver as DyAux uses *fastFVA* for efficiency.

---

## How to Use  
1. Clone the repository:  
   ```bash
   git clone https://github.com/RamanLab/DyAux.git
