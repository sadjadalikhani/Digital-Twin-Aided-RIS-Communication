function [W1,W2,epsilon1,epsilon2,flag] = ...
          wOptRobust(nRIS,nBS,nUE,G2,H3,Eprime,sigma,H1, ...
          rate,G1,h1_2path,emprCov)
rho = .05;
invT = load('./Channels/invT.mat').invT;
cvx_begin quiet
    variable W1(nBS,nBS) hermitian semidefinite
    variable W2(nBS,nBS) hermitian semidefinite
    variable epsilon1
    variable epsilon2
    variable x
    variable y
    expressions constr_1 constr_2 constr_3 intrfPow2 desiredPow1 desiredPow2 intrfPow2beamSweep
    intrfPow1 = 0;
    intrfPow2 = 0;
    desiredPow2 = 0; 
    
    U = invT'*kron(W1,1)*invT;
    u = invT*vec(h1_2path*W1);
    constant = trace(double(h1_2path'*h1_2path)*W1) - ...
        epsilon1*(2^(rate)-1);
    
    constr_1 = trace(U) - sqrt(2*log(1/rho))*x + ...
        log(rho)*y + constant;
    constr_2 = norm([vec(U);2*u],2) - x;
    constr_3 = y*eye(nBS) + U;
    
    for i = 1:nUE
        A1 = double([H1(i,:); diag(G1(i,:))*H3]);
        A2 = double([zeros(1,nBS); diag(G2(i,:))*H3]);
        intrfPow1 = intrfPow1 + real(trace(A2'*Eprime*A2*W1));
        intrfPow2 = intrfPow2 + real(trace(A1'*Eprime*A1*W2));
        desiredPow2 = desiredPow2 + real(trace(A2'*Eprime*A2*W2));
    end
    minimize(trace(W1)+trace(W2))
    subject to
        intrfPow1 <= epsilon2 - sigma^2;
        intrfPow2 <= epsilon1 - sigma^2;
        real(constr_1) >= 0;
        real(constr_2) <= 0;
        constr_3 == hermitian_semidefinite(nBS);
        y >= 0;
        desiredPow2 >= epsilon2*(2^rate-1);
        epsilon1 >= sigma^2;
        epsilon2 >= sigma^2;
cvx_end

% [W1] = rankOneApprox(W1,nBS,nRIS,nUE,G1,G2,H3,sigma,rate, ...
%                 epsilon1,epsilon2,W1,W2,H1,3,3);
% [W2] = rankOneApprox(W2,nBS,nRIS,nUE,G1,G2,H3,sigma,rate, ...
%                 epsilon1,epsilon2,W1,W2,H1,3,3);

% [W1] = gaussRand(W1)
% [W2] = gaussRand(W2)
            
if (isnan(W1) | isnan(W2) | isnan(epsilon1) | isnan(epsilon2))
    flag = 1;
else
    flag = 0;
end
end