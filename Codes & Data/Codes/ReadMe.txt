

**** Codes for the algorithm 
ADMM_Optim.R        --   Function on ADMM algorithm to estimate the parameters in our model given network data.
                    --   Take input as adjacency matrix of network data X, and tuning parameters gamma and delta.

GD.R                --   Gradient descent function used for updating X_m^alpha and X_m^M in step 1 of ADMM algorithm

**** Codes for generating data and running numerical experiments on synthetic data
SynData.R           --   Function for creating random graph (synthetic data) with the input parameters alpha, F, D, S, N
Synthetic Data.R    --   Codes to set for generating alpha, F, D, S, and for the procedure to choose tuning parameters.
                    --   to reproduce the results, please start with this file and you will be prompted run other code when 
                         necessary

**** Codes for the real data analysis on statistician's network data
citation_network_analysis.R   --  For reproducing the results in data analysis part of real statistician's data 
                                  presented in paper