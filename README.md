# COP
Optimization models to meet the conditions of order preservation in the analytic hierarchy process

# Installation
Before attempting to apply the MNVLLSM and NPRAOC, make sure your system has the following software available: 

[matlab 2018b or higher](https://www.mathworks.com/products/matlab.html)

[yalmip](https://yalmip.github.io/)

[Ipopt](https://github.com/coin-or/Ipopt)

[gurobi 10 or higher](https://www.gurobi.com/)


# NV_index_exchangeability

[NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(A), where A is the pairwise comparison matrix provided by the decision maker. This function
aims to calculate the number of violations and indentify the pairwise comparisons violating the index-exchangeability. 

The output parameter result of this function is as follows：

| parameter               |                          |
| --------------------- | ---------------------------- |
| NV                               | The total number of violations   |
| violation_index_equal            | The pairwise comparisons violate $$ a_{ij} = 1 \Leftrightarrow w_i=w_j$$ |
| violation_index_inequal          | 模型中的隐藏层数量           |

$$ a_{ij} = 1 \Leftrightarrow w_i=w_j$$
