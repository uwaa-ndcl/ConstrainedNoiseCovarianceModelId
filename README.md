# Autocovariance Least Squares with Constrained Noise Covariance Model Identification

---

The scripts to run the constrained Autocovariance Least Squares technique for two examples are contained in this repository.

## Mass Spring Damper Example
![image](https://github.com/user-attachments/assets/f315cbe3-5d15-4c0f-ab0a-f57451c6ab91)


1. [gen_data_msd](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MassSpringDamper/gen_data_msd.m)
   - Creates simulated datasets of the msd dynamics with process and measurement noise
   - Generates ALS inputs with initial suboptimal process noise covariance, $Q$ and measurement noise covariance,  $R$.
2. [run_als_msd](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MassSpringDamper/run_als_msd.m)
   - Run script for constrained ALS problem
   - Calls [setup_ALS_msd.m](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MassSpringDamper/setup_ALS_msd.m) with defines lags and other ALS inputs
   - [als_msd](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MassSpringDamper/als_msd.m)
     - ALS class with mass spring damper constraints for 7 temperature problem
3. [plot_lags_msd](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MassSpringDamper/plot_lags_msd.m)
   - Plots results from constrained ALS problem
   - Saves mean and standard deviation of $Q$ and $R$ solutions
   - Figure 2
4. [plot_QR_T](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MassSpringDamper/plot_QR_T.m)
   - Figure 5
5. [Phi_f](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MassSpringDamper/Phi_F.m): Figure 3
  
## Model for Aeroelastic Response to Gust Excitation
![image](https://github.com/user-attachments/assets/ad53d107-e6d3-4d9d-9e46-7aec70e59f63)

1. [wtData_setup](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MARGE/wtData_setup.m)
   - Loads wind tunnel datasets: 
   - Relies on models located here: 
   - Generates ALS inputs with initial suboptimal process noise covariance, $Q$ and measurement noise covariance,  $R$.
2. [run_als_MARGE](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MARGE/run_als_MARGE.m)
   - Run script for constrained ALS problem
   - Calls [setup_ALS_MARGE.m](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MARGE/setup_ALS_MARGE.m) with defines lags and other ALS inputs
   - [als_MARGE](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MARGE/als_MARGE.m)
     - ALS class with MARGE constraints for 5 dynamic pressure problem
3. [plot_lags_MARGE](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MARGE/plot_lags_MARGE.m)
   - Plots results from constrained ALS problem
   - Saves mean and standard deviation of $Q$ and $R$ solutions
4. [plot_QR_qbar](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MARGE/plot_QR_qbar.m)
   - Figure 6
   - Relies on unconstrained ALS solutions:d
5. [wtData_ALSest](https://github.com/uwaa-ndcl/ConstrainedNoiseCovarianceModelId/blob/main/MARGE/wtData_ALSest.m)
   - Figures 7-10
   - Relies on wind tunnel data: 

