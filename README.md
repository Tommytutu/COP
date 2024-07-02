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

  Apply the function **[NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(A)**, we get $NV=54$. Thus, $\mathbf{A}$ violates the index-exchangeability condition.
  In such cases, we have two options: directly derive the priority vector by the function **MNVLLSM**, goes to Step 2; communicate with the DM and use function **NPRAOC** to provide some
  modification suggestions and get more coherent preferences, goes to Step 3.

  ```matlab
  clear;
  clc;
  A=[1    	5    	3    	7    	6    	6    	 1/3	 1/4
   1/5	1    	 1/3	5    	3    	3    	 1/5	 1/7
   1/3	3    	1    	6    	3    	4    	6    	 1/5
   1/7	 1/5	 1/6	1    	 1/3	 1/4	 1/7	 1/8
   1/6	 1/3	 1/3	3    	1    	 1/2	 1/5	 1/6
   1/6	 1/3	 1/4	4    	2    	1    	 1/5	 1/6
   3    	5    	 1/6	7    	5    	5    	1    	 1/2
   4    	7    	5    	8    	6    	6    	2    	1    ];


  [NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(A);
  
  
  
* **Step 2**: Derive the priority vector by function **MNVLLSM**
  ```matlab
  clear;
  clc;
  A=[1    	5    	3    	7    	6    	6    	 1/3	 1/4
   1/5	1    	 1/3	5    	3    	3    	 1/5	 1/7
   1/3	3    	1    	6    	3    	4    	6    	 1/5
   1/7	 1/5	 1/6	1    	 1/3	 1/4	 1/7	 1/8
   1/6	 1/3	 1/3	3    	1    	 1/2	 1/5	 1/6
   1/6	 1/3	 1/4	4    	2    	1    	 1/5	 1/6
   3    	5    	 1/6	7    	5    	5    	1    	 1/2
   4    	7    	5    	8    	6    	6    	2    	1    ];
  [NV,GCI,weight,sotime]=MNVLLSM(A)

  results
  Optimize a model with 6497 rows, 2600 columns and 24600 nonzeros
  Model fingerprint: 0x6ee263a4
  Model has 36 quadratic objective terms
  Variable types: 8 continuous, 2592 integer (2592 binary)
  Coefficient statistics:
    Matrix range     [1e+00, 4e+01]
    Objective range  [2e-01, 4e+01]
    QObjective range [2e-01, 7e-01]
    Bounds range     [1e+00, 1e+00]
    RHS range        [1e-02, 4e+01]
 Found heuristic solution: objective 14435.126805
 Presolve removed 3931 rows and 1566 columns
 Presolve time: 0.02s
 Presolved: 2566 rows, 1034 columns, 11796 nonzeros
 Presolved model has 36 quadratic objective terms
 Variable types: 8 continuous, 1026 integer (1026 binary)

 Root relaxation: objective 3.903616e+00, 2012 iterations, 0.03 seconds (0.03 work units)

     Nodes    |    Current Node    |     Objective Bounds      |     Work
  Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

      0     0    3.90362    0  125 14435.1268    3.90362   100%     -    0s
 H    0     0                    1866.1973611    3.90362   100%     -    0s
      0     0  495.35402    0   72 1866.19736  495.35402  73.5%     -    0s
 H    0     0                    1210.1388818  495.76171  59.0%     -    0s
      0     0  495.76171    0   65 1210.13888  495.76171  59.0%     -    0s
      0     0  495.76171    0   66 1210.13888  495.76171  59.0%     -    0s
      0     0  536.39241    0   65 1210.13888  536.39241  55.7%     -    0s
      0     2  567.11167    0   65 1210.13888  567.11167  53.1%     -    0s
 H  126    63                    1169.1407201  592.03409  49.4%  58.9    0s
 H  191    84                    1087.9426931  592.03409  45.6%  47.7    0s
 H  274   106                    1087.1517747  619.57017  43.0%  54.9    0s
 H  619   246                    1066.5582293  620.04202  41.9%  42.7    0s
 H  621   246                    1025.5570591  620.04202  39.5%  42.7    0s
 H  959   335                    1005.2893516  741.12853  26.3%  37.4    1s
 H  960   335                    1005.2344244  741.12853  26.3%  37.4    1s
 H  961   335                    1005.2212685  741.12853  26.3%  37.4    1s
 H  962   335                    1005.2159270  741.12853  26.3%  37.3    1s
 H  963   335                    1005.1766118  741.12853  26.3%  37.3    1s
 H 1214   345                    1005.1114567  769.08480  23.5%  35.0    1s
 H 1215   345                    1005.0580807  769.08480  23.5%  34.9    1s

 Cutting planes:
   Learned: 82
   Gomory: 14
   Flow cover: 32
   RLT: 3

 Explored 4106 nodes (128686 simplex iterations) in 4.59 seconds (3.51 work units)
 Thread count was 12 (of 12 available processors)

 Solution count 9: 1005.06 1005.11 1025.56 ... 14435.1

 Optimal solution found (tolerance 1.00e-04)
 Best objective 1.005058080674e+03, best bound 1.005057404310e+03, gap 0.0001%

 NV =

    24.5000


 GCI =

     0.5581


 weight =

     0.1724    0.0576    0.1212    0.0193    0.0409    0.0409    0.1776    0.3699


 sotime =

     4.6278


  

  
