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

 

  %% The following is the solving process and results.

  % CPU model: Intel(R) Core(TM) i7-9750H CPU @ 2.60GHz, instruction set [SSE2|AVX|AVX2]
  % Thread count: 6 physical cores, 12 logical processors, using up to 12 threads
  
  % Optimize a model with 6497 rows, 2600 columns and 24600 nonzeros
  % Model fingerprint: 0x6ee263a4
  % Model has 36 quadratic objective terms
  % Variable types: 8 continuous, 2592 integer (2592 binary)
  % Coefficient statistics:
  %   Matrix range     [1e+00, 4e+01]
  %   Objective range  [2e-01, 4e+01]
  %   QObjective range [2e-01, 7e-01]
  %   Bounds range     [1e+00, 1e+00]
  %   RHS range        [1e-02, 4e+01]
  % Found heuristic solution: objective 14435.126805
  % Presolve removed 3931 rows and 1566 columns
  % Presolve time: 0.02s
  % Presolved: 2566 rows, 1034 columns, 11796 nonzeros
  % Presolved model has 36 quadratic objective terms
  % Variable types: 8 continuous, 1026 integer (1026 binary)
  
  % Root relaxation: objective 3.903616e+00, 2012 iterations, 0.03 seconds (0.03 work units)
  
  %     Nodes    |    Current Node    |     Objective Bounds      |     Work
  %  Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time
  % 
  %      0     0    3.90362    0  125 14435.1268    3.90362   100%     -    0s
  % H    0     0                    1866.1973611    3.90362   100%     -    0s
  %      0     0  495.35402    0   72 1866.19736  495.35402  73.5%     -    0s
  % H    0     0                    1210.1388818  495.76171  59.0%     -    0s
  %      0     0  495.76171    0   65 1210.13888  495.76171  59.0%     -    0s
  %      0     0  495.76171    0   66 1210.13888  495.76171  59.0%     -    0s
  %      0     0  536.39241    0   65 1210.13888  536.39241  55.7%     -    0s
  %      0     2  567.11167    0   65 1210.13888  567.11167  53.1%     -    0s
  % H  126    63                    1169.1407201  592.03409  49.4%  58.9    0s
  % H  191    84                    1087.9426931  592.03409  45.6%  47.7    0s
  % H  274   106                    1087.1517747  619.57017  43.0%  54.9    0s
  % H  619   246                    1066.5582293  620.04202  41.9%  42.7    0s
  % H  621   246                    1025.5570591  620.04202  39.5%  42.7    0s
  % H  959   335                    1005.2893516  741.12853  26.3%  37.4    1s
  % H  960   335                    1005.2344244  741.12853  26.3%  37.4    1s
  % H  961   335                    1005.2212685  741.12853  26.3%  37.4    1s
  % H  962   335                    1005.2159270  741.12853  26.3%  37.3    1s
  % H  963   335                    1005.1766118  741.12853  26.3%  37.3    1s
  % H 1214   345                    1005.1114567  769.08480  23.5%  35.0    1s
  % H 1215   345                    1005.0580807  769.08480  23.5%  34.9    1s
   
  % Cutting planes:
  %   Learned: 82
  %   Gomory: 14
  %   Flow cover: 32
  %   RLT: 3
   
  % Explored 4106 nodes (128686 simplex iterations) in 4.59 seconds (3.51 work units)
  % Thread count was 12 (of 12 available processors)
  
  % Solution count 9: 1005.06 1005.11 1025.56 ... 14435.1
   
  % Optimal solution found (tolerance 1.00e-04)
  % Best objective 1.005058080674e+03, best bound 1.005057404310e+03, gap 0.0001%
   
  % NV =
  
  %    24.5000
  
   
  % GCI =
   
  %     0.5581
  
  
  % weight =
  % 
  %     0.1724    0.0576    0.1212    0.0193    0.0409    0.0409    0.1776    0.3699
   
  
  % sotime =
   
  %     4.6278
  
* **Step 3**: use function **NPRAOC** to provide some
  modification suggestions and get more coherent preferences

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

  n=length(A);
  % determine the inconsistency index threshold
  if n==3
        CI_bar=0.31
  elseif n==4
        CI_bar=0.35
  else
        CI_bar=0.37
  end

  %determine the elements that the DM refuses to  revise if b(i,j)=1,then the elements a(i,j)should no be revised. For example, if the DM refuses to revise the element a(3,7), then we can set b(3,7)=1
  
  b=zeros(n,n)
  b(3,7)=1 

  [a_bar, weight]=NPRAOC(A, CI_bar,b)

  %% The following is the solving process and results.

  % Optimize a model with 5738 rows, 2697 columns and 23584 nonzeros
  % Model fingerprint: 0x868e0566
  % Model has 1 quadratic constraint
  % Variable types: 165 continuous, 2532 integer (2532 binary)
  % Coefficient statistics:
  %   Matrix range     [1e-06, 1e+04]
  %   QMatrix range    [1e+00, 1e+00]
  %   Objective range  [2e-02, 1e+02]
  %   Bounds range     [1e+00, 1e+00]
  %   RHS range        [6e-01, 1e+04]
  % Presolve removed 2344 rows and 942 columns
  % Presolve time: 0.19s
  % Presolved: 3395 rows, 1818 columns, 15888 nonzeros
  % Presolved model has 63 quadratic constraint(s)
  % Variable types: 187 continuous, 1631 integer (1631 binary)
 
  % Root relaxation: objective 1.043332e+02, 937 iterations, 0.03 seconds (0.05 work units)
  % 
  %     Nodes    |    Current Node    |     Objective Bounds      |     Work
  %  Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time
  % 
  %      0     0  104.33320    0  191          -  104.33320      -     -    0s
  %      0     0  204.33096    0  202          -  204.33096      -     -    0s
  %      0     0  204.33096    0  203          -  204.33096      -     -    0s
  %      0     0  603.90462    0  187          -  603.90462      -     -    0s
  %      0     0  603.91078    0  188          -  603.91078      -     -    0s
  %      0     0  780.89906    0  198          -  780.89906      -     -    0s
  %      0     0  780.89906    0  199          -  780.89906      -     -    0s
  %      0     0  805.85128    0  214          -  805.85128      -     -    0s
  %      0     0  805.85128    0  214          -  805.85128      -     -    0s
  %      0     0  900.21484    0  202          -  900.21484      -     -    0s
  %      0     0  900.21484    0  206          -  900.21484      -     -    0s
  %      0     0  900.22165    0   98          -  900.22165      -     -    0s
  %      0     0  900.23198    0  174          -  900.23198      -     -    1s
  %      0     0  900.29279    0  191          -  900.29279      -     -    1s
  %      0     0  901.53324    0  248          -  901.53324      -     -    1s
  %      0     0  919.05726    0  299          -  919.05726      -     -    1s
  %      0     0  919.05726    0  300          -  919.05726      -     -    1s
  %      0     0  925.20244    0  247          -  925.20244      -     -    1s
  %      0     0  931.77213    0  287          -  931.77213      -     -    1s
  %      0     0  931.77213    0  321          -  931.77213      -     -    1s
  %      0     0  931.77213    0  319          -  931.77213      -     -    1s
  %      0     0  931.77213    0  326          -  931.77213      -     -    1s
  %      0     0  931.77213    0  329          -  931.77213      -     -    1s
  %      0     0  931.77213    0  339          -  931.77213      -     -    1s
  %      0     0  931.77213    0  338          -  931.77213      -     -    1s
  %      0     0  931.77213    0  320          -  931.77213      -     -    1s
  %      0     0  932.18720    0  125          -  932.18720      -     -    2s
  %      0     0  932.18720    0  128          -  932.18720      -     -    2s
  %      0     0  932.18720    0  237          -  932.18720      -     -    2s
  %      0     0  932.18720    0  199          -  932.18720      -     -    2s
  %      0     0  950.27055    0  304          -  950.27055      -     -    2s
  %      0     0  950.28163    0  263          -  950.28163      -     -    2s
  %      0     0  950.31442    0  256          -  950.31442      -     -    2s
  %      0     0  950.37069    0  279          -  950.37069      -     -    2s
  %      0     0  967.66395    0  260          -  967.66395      -     -    3s
  %      0     0  975.27761    0  253          -  975.27761      -     -    3s
  %      0     0 1001.73922    0  366          - 1001.73922      -     -    3s
  %      0     0 1003.29018    0  334          - 1003.29018      -     -    3s
  %      0     0 1003.29018    0  334          - 1003.29018      -     -    3s
  %      0     0 1012.64038    0  273          - 1012.64038      -     -    3s
  %      0     0 1012.64038    0  272          - 1012.64038      -     -    3s
  %      0     0 1017.21840    0  302          - 1017.21840      -     -    3s
  %      0     0 1033.03387    0  294          - 1033.03387      -     -    3s
  %      0     0 1035.05598    0  290          - 1035.05598      -     -    3s
  %      0     0 1039.93727    0  316          - 1039.93727      -     -    3s
  %      0     0 1041.46344    0  320          - 1041.46344      -     -    4s
  %      0     0 1042.11738    0  299          - 1042.11738      -     -    4s
  %      0     0 1044.56803    0  347          - 1044.56803      -     -    4s
  %      0     0 1044.58309    0  335          - 1044.58309      -     -    4s
  %      0     2 1047.80282    0  300          - 1047.80282      -     -    4s
  %     69    67 1500.30152   12   79          - 1104.81985      -   174    5s
  %    667   324 1400.46911   14  245          - 1177.18746      -  62.5   10s
  %    689   339 2000.46026   35  346          - 1186.39587      -  60.5   15s
  % H  888   397                    1500.3725781 1195.45508  20.3%  97.7   19s
  %    931   397 1412.39269   20   69 1500.37258 1195.63769  20.3%   104   20s
  %   1381   326 1374.72908   19   97 1500.37258 1290.86849  14.0%   123   25s
  %   1875   244     cutoff   23      1500.37258 1320.28995  12.0%   141   30s
  %   2587   118     cutoff   27      1500.37258 1400.51125  6.66%   150   35s
  % 
  % Cutting planes:
  %   Learned: 23
  %   Gomory: 26
  %   Cover: 422
  %   Projected implied bound: 4
  %   Clique: 122
  %   MIR: 152
  %   Mixing: 4
  %   StrongCG: 75
  %   Flow cover: 450
  %   GUB cover: 150
  %   Inf proof: 8
  %   Zero half: 133
  %   RLT: 7
  %   Relax-and-lift: 42
  % 
  % Explored 2828 nodes (428417 simplex iterations) in 35.76 seconds (35.93 work units)
  % Thread count was 12 (of 12 available processors)
  % 
  % Solution count 1: 1500.37 
  % 
  % Optimal solution found (tolerance 1.00e-04)
  % Best objective 1.500372578095e+03, best bound 1.500372578095e+03, gap 0.0000%
  % 
  % a_bar =
  
  %     1.0000    5.0000    0.1429    6.0000    4.0000    3.0000    0.3333    0.2500
  %     0.2000    1.0000    0.1111    3.0000    0.5000    0.3333    0.1667    0.1429
  %     7.0000    9.0000    1.0000    9.0000    8.0000    7.0000    6.0000    5.0000
  %     0.1667    0.3333    0.1111    1.0000    0.3333    0.2500    0.1429    0.1250
  %     0.2500    2.0000    0.1250    3.0000    1.0000    0.5000    0.2000    0.1667
  %     0.3333    3.0000    0.1429    4.0000    2.0000    1.0000    0.2500    0.2000
  %     3.0000    6.0000    0.1667    7.0000    5.0000    4.0000    1.0000    0.5000
  %     4.0000    7.0000    0.2000    8.0000    6.0000    5.0000    2.0000    1.0000
 
  % weight =
  
  %     0.0697    0.0265    0.4818    0.0192    0.0366    0.0505    0.1327    0.1832

