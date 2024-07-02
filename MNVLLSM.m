
function [NV,GCI,weight,sotime]=MNVLLSM(a)
n=length(a);
M=41;
e=10^(-2);

la=[];
ea=[];
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                if a(i,j)>a(k,l)
                    la(i,j,k,l)=1;
                elseif a(i,j)==a(k,l)
                    ea(i,j,k,l)=1;
          
                end
            end
        end
    end
end



% decision variables
y=sdpvar(1,n);
lw=binvar(n,n,n,n,'full');
ew=binvar(n,n,n,n,'full');

% onjective functions
% stage 1, minimize NV
ob1=0;


if max(max(a))<=1
     ob1=0;
else
    for i=1:n
    for j=i+1:n
        for k=1:n
            for l=k+1:n
                ob1=ob1+la(i,j,k,l)*(1-lw(i,j,k,l))+0.25*ea(i,j,k,l)*(1-ew(i,j,k,l));
            end
        end
    end
    end
end

    



% Stage 2
ob2=0;

for i=1:n
    for j=1:n
        ob2=ob2+(log(a(i,j))-y(i)+y(j))^2/((n-2)*(n-1));
       
    end
end

object=M*ob1 + ob2;

co=[];
for i=1:n
    for j=i:n
        for k=1:n
            for l=k:n
                co=[co,y(i)+y(l)-(y(k)+y(j))+M*(1-lw(i,j,k,l))-e>=0];
                co=[co,y(i)+y(l)-(y(k)+y(j))-M*lw(i,j,k,l)<=0];
                co=[co,y(i)+y(l)-(y(k)+y(j))+M*(1-ew(i,j,k,l))>=0];
                co=[co,y(i)+y(l)-(y(k)+y(j))-M*(1-ew(i,j,k,l))<=0];
                co=[co,y(i)+y(l)-(y(k)+y(j))-M*(lw(i,j,k,l)+ew(i,j,k,l))+e<=0];
            end
        end
    end
end




co=[co,sum(y)==0];

co=[co,-10<=y<=10];
ops = sdpsettings('solver', 'gurobi', 'gurobi.TuneTimeLimit', 0);



sol=optimize(co,object,ops);
sotime=sol.solvertime;

object=double(object);
ob1=double(ob1);
ob2=double(ob2);
NV=ob1;
GCI=ob2;
y=double(y);
lw=double(lw);
ew=double(ew);


w=exp(y);
weight=w./sum(w);  
for i=1:n
    for j=1:n
        W(i,j)=w(i)/w(j);
    end
end
end

