<h1 align="center"> Latent Factor + Sparse Matrix logistic regression model for network data </h1>

<p align="center">
  <a href="#overview">Overview</a> •
  <a href="#quickstart-with-the-data--models">Quickstart</a> •
  <a href="#acknowledgements">Acknowledgements</a> 
</p>

# Overview

- **Model Description** : Our paper proposes a combined model, which integrates the latent factor model and the logistic regression model, for the network data. It is noticed that neither a latent factor model nor a logistic regression model alone is sufficient to capture the structure of the data. The proposed model has a latent (i.e., factor analysis) model to represent the main technological trends (a.k.a., factors), and adds a sparse component that captures the remaining ad hoc dependence.

- **[Paper link](https://arxiv.org/abs/1912.00524)**: "Latent Factor + Sparse Matrix logistic regression model for network data"

- **Main functions**
    1. [ADMM_Optim.R](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/ADMM_Optim.R) : Main function for making inference for model parameters. Takes adjacency matrix and tuning parameters of models (gamma,delta) as input parameters. Gives the estimated alpha, matrices M, L and S as output of the function.
    2. [Model_Sel.R](https://github.com/namjoonsuh/Citation-Network/blob/master/Codes%20%26%20Data/Codes/Synthetic%20Networks.R) : Given the adjacency matrix of the network data and the ranges of grids for gamma and delta to search over, the function gives: 
      **(1)** a pair of indices (gamma, delta) that minimizes AIC over the given grid. 
      **(2)** a pair of indices (gamma, delta) that minimizes BIC over the given grid. 
      **(3)** the number of non-zero entries of the estimated S for each point on the grid. 
      **(4)** the rank of the estimated L for each point on the grid.
    
# Karate club data example

```R
import load_data
# first time it runs, downloads and caches the data
df = load_data.load_county_level(data_dir='/path/to/data') 
```
