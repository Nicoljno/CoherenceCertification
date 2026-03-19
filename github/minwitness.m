function rho = minwitness(d,m,W,channel,cc)

%-------------------------------------------------------------------------%
%This function evaluates the states which minimise the witness after
%passing through a channel \Lambda.
%
%Input:
%- d: dimension of the states;
%- m: number of states in the set;
%- W: witness operators stored as cells W{x} for x = 1,...,m;
%- channel: '1' for depolarising channel, '2' for amplitude-damping
%channel;
%- cc: coefficients for amplitude-damping channel stored as cells c{j,i}
%for i = 1,...,d and j = i+1,...,d;
%
%Output:
%- rho: states which minimise the witness.
%-------------------------------------------------------------------------%

%Identity
id = eye(d);

if channel == 1
    %Witness operators
    for j = 1 : m
        [V,D] = eig(W{j});
        [maxlambda,ind] = min(diag(D));
        rho{j} = V(:,ind)*V(:,ind)';
    end
elseif channel == 2
    %Maximally entangled state
    phi = MaxEntangled(d);

    %Amplitude-damping channel
    K = ampdamping(d,cc);

    %Choi state
    eta = 0;
    for j = 1 : d*(d-1)/2+1
        eta = eta + Tensor(K{j},id)*phi*phi'*Tensor(K{j}',id);
    end

    %Witness operators
    for x = 1 : m
        beta{x} = transpose(d*PartialTrace(Tensor(W{x},id)*eta,1,[d d]));
        [V,D] = eig(beta{x});
        [minlambda,ind] = min(diag(D));
        rho{x} = V(:,ind)*V(:,ind)';
    end
end

end