# COP
Optimization models to meet the conditions of order preservation in the analytic hierarchy process

# Installation
Before attempting to apply the MNVLLSM and NPRAOC, make sure your system has the following software available: 

* [matlab 2018b or higher](https://www.mathworks.com/products/matlab.html)

* [yalmip](https://yalmip.github.io/)

* [Ipopt](https://github.com/coin-or/Ipopt)

* [gurobi 10 or higher](https://www.gurobi.com/)


# NV_index_exchangeability

*[NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(A)*, where A is the pairwise comparison matrix provided by the decision maker. This function
aims to calculate the number of violations and indentify the pairwise comparisons violating the index-exchangeability. 

The output parameter result of this function is as follows：

* *NV*: the total number of violations;

* *violation_index_equal*: the pairwise comparisons violate $a_{ij}= a_{kl}  \Leftrightarrow a_{ik}= a_{jl}$

* *violation_index_inequal*: the pairwise comparisons violate $a_{ij}> a_{kl}  \Leftrightarrow a_{ik}> a_{jl}$


# MNVLLSM

*MNVLLSM* is a two-stage optimization model to derive the PCM's priority vector, where the first stage aims to minimize the numberof violation and the second stage minimizes the deviations  measured by the Logarithmic Least Squares Method (LLSM).

*[NV,GCI,weight,sotime]=MNVLLSM(A)*, the output parameter result of this function is as follows：

* *NV*: the total number of violations;
* *GCI*: the deviations measured bpy the LLSM;
* *weight*: the priority vector of PCM A with the minimal number of violations;
* *sotime*: the solving time of this model

# NPRAOC

*NPRAOC* aims to meet  index-exchangeability condition by revising some preferences.

*[**A_bar**, weight]=NPRAOC(**A**, CI_bar,**B**)*, the input and out parameters of this function is as follows：
* **A**: the PCM provided by the DM, $\mathbf{A}=(a_{ij})_{n\times n}$;
* *CI_bar*: the predetermined consistency index threshold for the Geometric Consistency Index (GCI), *CI_bar=0.31* when $n=3$; *CI_bar=0.35* when $n=4$ and *CI_bar=0.37* when $n>4$;
* ***B***:  $\mathbf{B}=(b_{ij})_{n\times n}$; if  the DM the DM refuses to modify the preference $$a_{ij}$$, then $b_{ij}=1$; otherwise, $b_{ij}=0$;
* ***A_bar***: the modified PCM that meets both the   index-exchangeability condition
and the acceptable inconsistency level;
* *weight*: the priority vector of the modified PCM ***A_bar*** that meets the POIP conditions.


# Example



