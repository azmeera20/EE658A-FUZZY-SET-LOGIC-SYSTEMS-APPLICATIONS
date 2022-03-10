%% Divide fuzzy regions into 2N+1 regions
function [ out_Reg, N_R ] = FuzzyRegions( N, x )
N_R = 2*N + 1;
out_Reg = cell(1,1);
a = zeros(1,1);
b = zeros(1,1);
c = zeros(1,1);
temp1 = (min(x));
temp2 = (max(x));
temp12 = (temp2 - temp1)/N;
div = temp12/2;
    for i = 1:N_R
        if i==1
            a(i,1) = min(x);
            b(i,1) = min(x);
            c(i,1) = min(x) + div;
        elseif i==N_R
            a(i,1) = max(x) - div;
            b(i,1) = max(x);
            c(i,1) = max(x);
        else
            a(i,1) = min(x) + (i-2)*div;
            b(i,1) = min(x) + (i-1)*div;
            c(i,1) = min(x) + (i-0)*div;
        end
        
        if i==1
            out_Reg{i} = trapmf(x,[0,0,b(i,1),c(i,1)]);
        elseif i==N_R
            out_Reg{i} = trapmf(x,[a(i,1),b(i,1),c(i,1)+div,c(i,1)+div]);
        else
            out_Reg{i} = trimf(x,[a(i,1),b(i,1),c(i,1)]);
        end
    [temp_val, temp_indx] = max(out_Reg{i});
    out_Reg{i}(temp_indx) = ceil(temp_val);
    end
end