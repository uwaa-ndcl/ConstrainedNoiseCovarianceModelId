# Autocovariance Least Squares with Constrained Noise Covariance Model Identification

---

The scripts to run the constrained Autocovariance Least Squares technique for two examples are contained in this repository.

## Mass Spring Damper Example
![msdfbd]()

1. [gen_data_msd]()
   - Creates simulated datasets of the msd dynamics with process and measurement noise
   - Generates ALS inputs with initial suboptimal process noise covariance, $Q$ and measurement noise covariance,  $R$.
2. [run_als_msd]()
   - Run script for constrained ALS problem
   - Calls [setup_ALS_msd.m]() with defines lags and other ALS inputs
   - [als_msd]()
     - ALS class with mass spring damper constraints for 7 temperature problem
3. [plot_lags_msd]()
   - Plots results from constrained ALS problem
   - Saves mean and standard deviation of $Q$ and $R$ solutions
   - Figure 2
4. [plot_QR_T]()
   - Figure 5
5. [Phi_f.m](): Figure 3
  
## Model for Aeroelastic Response to Gust Excitation
![MARGE]()


1. [wtData_setup]()
   - Loads wind tunnel datasets: 
   - Relies on models located here: 
   - Generates ALS inputs with initial suboptimal process noise covariance, $Q$ and measurement noise covariance,  $R$.
2. [run_als_MARGE]()
   - Run script for constrained ALS problem
   - Calls [setup_ALS_MARGE.m]() with defines lags and other ALS inputs
   - [als_MARGE]()
     - ALS class with MARGE constraints for 5 dynamic pressure problem
3. [plot_lags_MARGE]()
   - Plots results from constrained ALS problem
   - Saves mean and standard deviation of $Q$ and $R$ solutions
4. [plot_QR_qbar]()
   - Figure 6
   - Relies on unconstrained ALS solutions:d
5. [wtData_ALSest]()
   - Figures 7-10
   - Relies on wind tunnel data: 

