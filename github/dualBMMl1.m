function [v,W] = dualBMMl1(d,m,rho,channel,cc)

%-------------------------------------------------------------------------%
%This function evaluates the dual of the BMM hierarchy at level 1 to
%certify the coherence of depolarising and amplitude-damping channels.
%
%Input:
%- d: dimension of the states;
%- m: number of states in the set;
%- rho: set of states stored as cells rho{x} for x = 1,...,m;
%- channel: '1' for depolarising channel, '2' for amplitude-damping
%channel;
%- cc: coefficients for amplitude-damping channel stored as cells c{j,i}
%for i = 1,...,d and j = i+1,...,d;
%
%Output:
%- v: critical visibility;
%- W: witness operators stored as cells W{x} for x = 1,...,m.
%-------------------------------------------------------------------------%

cvx_begin sdp quiet
cvx_solver mosek
variable Z(d*(m+1),d*(m+1)) hermitian semidefinite
variable ggamma(d,d,m+1,m+1) hermitian semidefinite
variable theta(d,d,m+1,m+1) hermitian semidefinite

%Channel definition
if channel == 1

elseif channel == 2
    %Amplitude-damping channel
    K = ampdamping(d,cc);

    for x = 1 : m
        sigma = 0;
        for j = 1 : d*(d-1)/2+1
            sigma = sigma + K{j}*rho{x}*K{j}';
        end
        rho{x} = sigma;
    end
end

%Depolarising channel
id = eye(d);
N = id/d;

%Indices of blocks
ind{1} = [1:d];
for k = 2 : m+1
    ind{k} = (k-1)*d + ind{1};
end

for i = 2 : m
    ZZ{i} = Z(ind{1},ind{i}) + Z(ind{i},ind{1}) + Z(ind{i},ind{i});
    gg{i} = 0;
    for j = i+1 : m+1
        gg{i} = gg{i} + ggamma(:,:,i,j);
    end
end
for j = 3 : m+1
    tt{j} = 0;
    for i = 2 : j-1
        tt{j} = tt{j} + theta(:,:,i,j);
    end
end
ZZ{m+1} = Z(ind{1},ind{m+1}) + Z(ind{m+1},ind{1}) + Z(ind{m+1},ind{m+1});

con = 1 + trace((ZZ{2}+gg{2})*(rho{1}-N)) + trace((ZZ{m+1}+tt{m+1})*(rho{m}-N));
for j = 3 : m
    con = con + trace((ZZ{j} + gg{j} + tt{j})*(rho{j-1}-N));
end
con == 0;

%Objective
obj = 1 + trace(Z(ind{1},ind{1})) + trace((ZZ{2}+gg{2})*rho{1}) + trace((ZZ{m+1}+tt{m+1})*rho{m});
for j = 3 : m
    obj = obj + trace((ZZ{j} + gg{j} + tt{j})*rho{j-1});
end

%Constraints
for i = 2 : m
    for j = i+1 : m+1
        Z(ind{i},ind{j}) + Z(ind{j},ind{i}) - ggamma(:,:,i,j) - theta(:,:,i,j) <= 0;
    end
end

minimize(real(obj))
cvx_end
v = cvx_optval;

%Witness operators
W{1} = ZZ{2} + gg{2};
for j = 2 : m-1
    W{j} = ZZ{j+1} + gg{j+1} + tt{j+1};
end
W{m} = ZZ{m+1} + tt{m+1};

end