function [a_bar, weight]=NPRAOC(a, CI_bar,b)

M=100;
e=0.001;
alph=100;

lisan=[ 1/9	 1/8	 1/7	 1/6	 1/5	 1/4	 1/3	 1/2	1    	2    	3    	4    	5    	6    	7    	8    	9];
duilisan=log(lisan)';

a=log(a);
n=length(a);

lx=binvar(n,n,n,n,'full');
x=sdpvar(n,n,'full');
derta=binvar(n,n,'full');

lijm=binvar(n,n,17,'full');
y=sdpvar(1,n);

%objective
distance=0;
for i=1:n
    for j=i+1:n
        distance=distance+abs((a(i,j))-x(i,j))/((n-1)*(n-2));
    end
end

NPR=0;
for i=1:n
    for j=i+1:n
       NPR=NPR+ derta(i,j);
    end
end


object=distance+alph*NPR;

%constraints
co=[];
for i=1:n
    for j=i:n
        for k=1:n
            for l=k:n
                co=[co,x(i,j)-x(k,l)+M*(1-lx(i,j,k,l))-e>=0];
                co=[co,x(i,j)-x(k,l)-M*lx(i,j,k,l)<=0];
                co=[co,(y(i)-y(j))-(y(k)-y(l))+M*(1-lx(i,j,k,l))-e>=0];
                co=[co,(y(i)-y(j))-(y(k)-y(l))-M*lx(i,j,k,l)<=0];
            end
        end
    end
end


    
for i=1:n
    for j=1:n
        co=[co,x(i,j)-a(i,j)<=M*derta(i,j)];
        co=[co,-x(i,j)+a(i,j)<=M*derta(i,j)];
    end
end


for i=1:n
    for j=1:n
        co=[co,x(i,j)+x(j,i)==0];
    end
end

for i=1:n
    for j=i+1:n
        co=[co,x(i,j)~=0];
    end
end

co=[co,sum(lijm,3)==1];
for i=1:n
    for j=1:n
       co=[co,x(i,j)-reshape(lijm(i,j,:),1,17)*duilisan==0];
    end
end

  ob2=0;
% 
for i=1:n
    for j=1:n
        ob2=ob2+(x(i,j)-y(i)+y(j)).^2/((n-2)*(n-1));
       
    end
end
co=[co,ob2<=CI_bar];
%co=[co,x(3,7)==log(6)];
for i=1:n
    for j=i:n
        if b(i,j)==1
            co=[co,x(i,j)==a(i,j)]
        end
    end
end

ops = sdpsettings('solver', 'gurobi', 'gurobi.TuneTimeLimit', 0);
sol=optimize(co,object,ops);
object=value(object);
ob2=value(ob2);
derta=value(derta);
NPR=value(NPR);
distance=value(distance);
x=value(x);
a_bar=exp(x);  
y=value(exp(y));
weight=y./sum(y);

end