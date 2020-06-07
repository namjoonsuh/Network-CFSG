<h1 align="center"> Latent Factor + Sparse Matrix logistic regression model for network data </h1>

<p align="center">
  <a href="#overview">Overview</a> •
  <a href="#quickstart-with-the-data--models">Quickstart</a> •
  <a href="#acknowledgements">Acknowledgements</a> 
</p>

# Overview

- ** Model Description ** : We propose a combined model, which integrates the latent factor model and the logistic regression model, for the network data.
It is noticed that neither a latent factor model nor a logistic regression model alone is sufficient to capture the structure of the data.
The proposed model has a latent (i.e., factor analysis) model to represents the main technological trends (a.k.a., factors), and adds a sparse component that captures the remaining ad-hoc dependence.

- **[Paper link](https://arxiv.org/abs/1912.00524)**: "Latent Factor + Sparse Matrix logistic regression model for network data"

# Karate club data example

```R
import load_data
# first time it runs, downloads and caches the data
df = load_data.load_county_level(data_dir='/path/to/data') 
```
