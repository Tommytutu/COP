function [NV,violation_index_equal, violation_index_inequal]=NV_index_exchangeability(a)

n=length(a)
la=[];
ea=[];
NV=0
violation_index_equal={}
violation_index_inequal={}
for i=1:n
    for j=i:n
        for k=1:n
            for l=k:n
                if a(i,j)>a(k,l) & a(i,k)<a(j,l)
                    NV=NV+1
                    violation_index_inequal=[violation_index_inequal, [i,j,k,l]]
                    
                elseif a(i,j)>a(k,l) & a(i,k)==a(j,l)
                    NV=NV+0.5
                    violation_index_equal=[violation_index_equal, [i,j,k,l]]
                else
                    NV=NV+0
                    
                end
             end
        end
    end
end
end