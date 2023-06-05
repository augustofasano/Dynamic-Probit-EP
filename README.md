# Expectation propagation for the smoothing distribution in dynamic probit

This repository is associated with the article [Anceschi, Fasano, Franzolini, and Rebaudo (2023)]() _Expectation propagation for the smoothing distribution in dynamic probit_. 
The key contribution of the paper is outlined below.

> we develop an expectation propagation (EP) approximation of the joint smoothing distribution in a dynamic probit model.

This repository provides codes to implement the inference methods associated with such a result. More precisely, we provide the `R` code to perform inference on the **smoothing distribution** under dynamic probit models with four different methods:

1. **i.i.d. sampling from the exact unified skew-normal distribution**
2. **expectation propagation (EP) approximation**
3. **partially factorized mean-field variational Bayes (PFM-VB) approximation**
4. **mean-field variational Bayes (MF-VB) approximation**

Structure of the repository:

* the functions to implement the above methods can be found in the `R` source file [`Functions.R`](https://github.com/augustofasano/Dynamic-Probit-EP/blob/main/Functions.R)
* the financial dataset analyzed in the illustration can be found in [`Financial.Rdata`](https://github.com/augustofasano/Dynamic-Probit-EP/blob/main/Financial.RData), while the original entire dataset is publicly available at [Yahoo Finance](https://finance.yahoo.com/)
* the code to reproduce the results in the paper is available at [`DynamicEP.R`](https://github.com/augustofasano/Dynamic-Probit-EP/blob/main/DynamicEP.R)