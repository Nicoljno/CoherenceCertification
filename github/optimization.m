clear all
clc

%-------------------------------------------------------------------------%
%Optimization algorithm to search for upper bounds on the critical
%visibility for a channel to be coherence-breaking.
%-------------------------------------------------------------------------%

%Dimension and number of states
d = 2;
m = 2;

%Identity
id = eye(d);

%Precision
epsilon = 10^-6;

%Channel

%Isotropic channel
%channel = 1; 

%Amplitude-damping channel
channel = 2;
for i = 1 : d
    for j = i+1 : d
        cc{j,i} = 0.5;
    end
end

%% Algorithm

%Starting ensemble
for x = 1 : m
    rho{x} = RandomDensityMatrix(d);
end

v0 = 0;
v = 1;
res = [];
while abs(v-v0) > epsilon
    %Iteration rule
    v0 = v;

    %Visibility
    [v,W] = dualBMMl1(d,m,rho,channel,cc);
    res = [res, v]

    %Maximum violation set
    rho = minwitness(d,m,W,channel,cc);

    %Save results
    save(['res_' num2str(d)],'res')

end