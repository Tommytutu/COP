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
Consider the Example 4 in the paper, $\mathbf{A}$ is a PCM with eight alternatives $\{x_1,x_2,x_3,x_4, x_5,x_6,x_7,x_8\}$ which has been discussed by  Saaty (2003).

$$\mathbf{A} =\left(\begin{array}{ccccccccccc}
 1     & 5     & 3     & 7     & 6     & 6     &  1/3  &  1/4 \\
     1/5  & 1     &  1/3  & 5     & 3     & 3     &  1/5  &  1/7 \\
     1/3  & 3     & 1     & 6     & 3     & 4     & 6     &  1/5 \\
     1/7  &  1/5  &  1/6  & 1     &  1/3  &  1/4  &  1/7  &  1/8 \\
     1/6  &  1/3  &  1/3  & 3     & 1     &  1/2  &  1/5  &  1/6 \\
     1/6  &  1/3  &  1/4  & 4     & 2     & 1     &  1/5  &  1/6 \\
    3     & 5     &  1/6  & 7     & 5     & 5     & 1     &  1/2 \\
    4     & 7     & 5     & 8     & 6     & 6     & 2     & 1     \\
  \end{array}\right)$$

* **Step 1**: Detect the violations of index-exchangeability condition

  Apply the function **[NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(A)**, we get $NV=5$. Thus, $\mathbf{A}$ violates the index-exchangeability condition.
  In such cases, we have two options: directly derive the priority vector by the function **MNVLLSM**, goes to Step 2; communicate with the DM and use function **NPRAOC** to provide some
  modification suggestions and get more coherent preferences, goes to Step 3.

  ```matlab
  clear;
  clc;
  A=[1    	5    	3    	7    	6    	4    	 1/3	 1/4
  1/5	1    	 1/3	4    	3    	 1/2	 1/6	 1/7
  1/3	3    	1    	6    	5    	3    	 1/4	 1/5
  1/7	 1/4	 1/6	1    	 1/2	 1/4	 1/8	 1/9
  1/6	 1/3	 1/5	2    	1    	 1/3	 1/7	 1/8
  1/4	2    	 1/3	4    	3    	1    	 1/5	 1/6
  3    	6    	4    	8    	7    	5    	1    	 1/2
  4    	7    	5    	9    	8    	6    	2    	1    ];

  [NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(A);


  
  The results are as follows:
  NV=5;

  violation_index_equal=([2,3,1,8],	[2,6,1,7],	[3,3,2,6],	[3,5,1,6],	[4,4,2,6],	[4,5,1,7],	[4,6,3,8],	[5,5,2,6],	[5,6,1,8],	[5,6,3,7]);
  violation_index_inequal={};


* **Step 2**: Derive the priority vector by function **MNVLLSM**
  

  
