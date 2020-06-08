<h1 align="center"> Latent Factor + Sparse Matrix logistic regression model for network data </h1>

<p align="center">
  <a href="#overview">Overview</a> •
  <a href="#quickstart-with-the-data--models">Quickstart</a> •
  <a href="#acknowledgements">Acknowledgements</a> 
</p>

# Overview

- **Model Description** : Our paper proposes a combined model, which integrates the latent factor model and the logistic regression model, for the network data. It is noticed that neither a latent factor model nor a logistic regression model alone is sufficient to capture the structure of the data. The proposed model has a latent (i.e., factor analysis) model to represent the main technological trends (a.k.a., factors), and adds a sparse component that captures the remaining ad hoc dependence.

- **[Paper link](https://arxiv.org/abs/1912.00524)**: "Latent Factor + Sparse Matrix logistic regression model for network data"

- **Main functions for model selection**
    1. [ADMM_Optim](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/ADMM_Optim.R) : Main function for making inference for model parameters. Takes adjacency matrix and tuning parameters of models (gamma,delta) as input parameters. Gives the estimated alpha, matrices M, L and S as output of the function.
    2. [Model_Sel](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Synthetic%20Networks.R) : Given the adjacency matrix of the network data and the ranges of grids for gamma and delta to search over, the function gives: 
      **(1)** a pair of indices (gamma, delta) that minimizes AIC over the given grid. 
      **(2)** a pair of indices (gamma, delta) that minimizes BIC over the given grid. 
      **(3)** the number of non-zero entries of the estimated S for each point on the grid. 
      **(4)** the rank of the estimated L for each point on the grid.
    3. [CV](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Synthetic%20Networks.R) : Code for Network Cross-validation. Refer Section 6.2. of the [paper](https://arxiv.org/abs/1912.00524) for detailed explanation of the procedure. Takes the adjacency matrix of the network data, a pair of tuning parameters (gamma, delta) and the number of K iterations for the averaged mis-classifation rate of edges. Gives the K averaged mis-classification rate as output. 
    4. [Eval_func](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Synthetic%20Networks.R) : Code for evaluation of the selected model. Given the adjacency matrix of the network data, and selected model parameters (gamma, delta) through AIC, BIC or Heuristic Network Cross-validation, the function gives:
    **(1)** the rank of estimated L matrix. (i.e. K)
    **(2)** K clustered nodes by applying k-means algorithm on K eigen-vectors of estimated L matrix. 
    **(3)** the number of non-zero entries on upper-triangular part of the estimated S matrix. 
    **(4)** a list of pairs of nodes which create the ad-hoc edges of the selected model. 
 
- **Codes for Simulations and Real Data analysis**
We made a simple demonstration on how our model can be applied to a network data by using aforementioned functions.
A famous [Zachary Karate Club Network](https://en.wikipedia.org/wiki/Zachary%27s_karate_club) data is analzyed [here](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Zhacary_Karate.R). 

    1. [Simulation1](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Simulation1.R) : Code for reproducing results presented in Section 6.1. of the paper.
    2. [Simulation2](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Simulation2.R) : Code for reproducing results presented in Section 6.2. of the paper.
    3. [Polblog](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Polblog.R) : Code for reproducing results presented in Section 7.1. of the paper.
    4. [Citation_network](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/citation_network.R) : Code for reproducing results presented in Section 7.2. of the paper.

- **Data References**




