# COP
Optimization models to meet the conditions of order preservation in the analytic hierarchy process

# Installation
Before attempting to apply the MNVLLSM and NPRAOC, make sure your system has the following software available: 

[matlab 2018b or higher](https://www.mathworks.com/products/matlab.html)

[yalmip](https://yalmip.github.io/)

[Ipopt](https://github.com/coin-or/Ipopt)

[gurobi 10 or higher](https://www.gurobi.com/)


# NV_index_exchangeability
<span style="text-shadow: 2px 2px 2px grey;">[NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(A)</span>
[NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(A), where A is the pairwise comparison matrix provided by the decision maker. This function
aims to calculate the number of violations and indentify the pairwise comparisons violating the index-exchangeability. 

The output parameter result of this function is as follows：

| parameter               |                          |
| --------------------- | ---------------------------- |
| NV         | The total number of violations   |
| Batch Size            | 每次迭代时用于训练的样本数量 |
| Number of Layers      | 模型中的隐藏层数量           |
| Hidden Units          | 每个隐藏层中的神经元数量     |
| Activation Function   | 用于引入非线性变换的函数     |
| Dropout Rate          | 在训练过程中随机丢弃部分神经元以防止过拟合 |
| Regularization Parameter | 用于控制模型复杂度的正`

 $$E=mc^2$$
