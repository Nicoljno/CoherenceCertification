% Block-moment-matrix relaxation for obtaining upper bounds on the critical
% visibility for equiangular tight frames (ETF) of dimension dim.
%
% The ETFs for dimension [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13] 
% are saved in the file "ETF_2d.mat"

clear all
clc

reslist=[];
dims=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];

for n_=2:length(dims)

d=dims(n_); N=2*d;

S = load("ETF_2d.mat");
name = sprintf('vec_%d', d);
S = S.(name);

id=eye(d);

A=zeros(d,d,N);

for k=1:N
    A(:,:,k)=S(:,k)*S(:,k)';
end
tick=0;
for k=1:N
    for m=1:k-1
        tick=tick+1;
        map(m,k)=tick;
    end
end
tic
nHerm = nchoosek(N,2);
cvx_begin sdp  quiet
cvx_solver mosek
variable v
variable M(d,d,nHerm) hermitian semidefinite


for k=1:N
    r{k}=v*A(:,:,k)+(1-v)/d*eye(d);
end


row0=[id];
for  k=1:N
    row0=[row0 r{k}];
end
G=[row0];

for m=1:N
    dum=[r{m}];
    for  k=1:N
        if k==m        
        dum=[dum r{k}];
        elseif m<k
            ind=map(m,k);
        dum=[dum M(:,:,ind)];    
        M(:,:,ind) <= r{m};
        M(:,:,ind) <= r{k};
        elseif m>k
            ind=map(k,m);
        dum=[dum M(:,:,ind)];
        end
    end
    G=[G;dum];
end
G>=0;
    
    

maximise(v)
cvx_end
vBMM=cvx_optval

time=toc

reslist=[reslist; vBMM time]
end
save('result_ETF','reslist')


%% DUAL formulation of the problem

cvx_begin sdp quiet
    cvx_solver mosek

    variable X(d*(N+1), d*(N+1)) hermitian semidefinite

    expression obj
    expression expr

    blk0 = 1:d;

    obj = trace( X(blk0,blk0) );
    for k = 1:N
        blkk = k*d + (1:d);
        obj = obj + (1/d) * real(trace( X(blk0,blkk) + X(blkk,blk0) + X(blkk,blkk) ));
    end
    minimize( obj );

    expr = 0;
    for k = 1:N
        blkk = k*d + (1:d);
        Bk = A(:,:,k) - eye(d)/d;
        expr = expr + real(trace( ( X(blk0,blkk) + X(blkk,blk0) + X(blkk,blkk) ) * Bk ));
    end
    expr == -1;

    for m = 1:N
        blkm = m*d + (1:d);
        for k = m+1:N
            blkk = k*d + (1:d);

            X(blkm,blkk) + X(blkk,blkm) <= 0;
        end
    end
cvx_end
cvx_optval

