clear;
clc;
a=[1    	5    	3    	7    	6    	4    	 1/3	 1/4
 1/5	1    	 1/3	4    	3    	 1/2	 1/6	 1/7
 1/3	3    	1    	6    	5    	3    	 1/4	 1/5
 1/7	 1/4	 1/6	1    	 1/2	 1/4	 1/8	 1/9
 1/6	 1/3	 1/5	2    	1    	 1/3	 1/7	 1/8
 1/4	2    	 1/3	4    	3    	1    	 1/5	 1/6
3    	6    	4    	8    	7    	5    	1    	 1/2
4    	7    	5    	9    	8    	6    	2    	1    ];

n=length(a);

% determine the inconsistency index threshold
if n==3
    CI_bar=0.31
elseif n==4
    CI_bar=0.35
else
    CI_bar=0.37
end

%%determine the elements that the DM refuses to  revise
%if b(i,j)=1,then the elements a(i,j)should no be revised
%For example, if the DM refuses to revise the element a(3,7), then we can
%set b(3,7)=1 
b=zeros(n,n)
b(3,7)=1 

[a_bar, weight]=NPRAOC(a,CI_bar,b)

% if the DM refuses to break his/her preference cycle, then NPRAOC will
% infeasible. Take the example of Example 4, if the DM insists that the
% preference cycle between alternatives 1,3,7 are reasonable, that is,
% b(1,7)=1, b(1,3)=1, b(3,7)=1, then NPRAOC will be infeasible.
b=zeros(n,n)
b(1,7)=1
b(1,3)=1
b(3,7)=1
[a_bar, weight]=NPRAOC(a,CI_bar,b)


[NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(a)

