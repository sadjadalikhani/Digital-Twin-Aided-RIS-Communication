function [W1,W2,epsilon1,epsilon2,flag] = ...
    wOpt(nRIS,nBS,nUE,G2,H3,Eprime,sigma,H1,rate,G1)
cvx_begin quiet
    variable W1(nBS,nBS) hermitian semidefinite
    variable W2(nBS,nBS) hermitian semidefinite
    variable epsilon1
    variable epsilon2
    expressions intrfPow1 intrfPow2 desiredPow1 desiredPow2 intrfPow2beamSweep
    intrfPow1 = 0;
    intrfPow2 = 0;
    desiredPow2 = 0;
    intrfPow2beamSweep = 0;
    
    desiredPow1 = real(trace(H1'*H1*W1));
    for i = 1:nUE
        A1 = double([H1(i,:); diag(G1(i,:))*H3]);
        A2 = double([zeros(1,nBS); diag(G2(i,:))*H3]);
        A3 = double([H1(i,:); zeros(nRIS,nBS)]);
        intrfPow1 = intrfPow1 + real(trace(A2'*Eprime*A2*W1));
        intrfPow2 = intrfPow2 + real(trace(A1'*Eprime*A1*W2));
        desiredPow2 = desiredPow2 + real(trace(A2'*Eprime*A2*W2));
        intrfPow2beamSweep = intrfPow2beamSweep + real(trace(A3*W2*A3'*Eprime));
    end
    minimize(trace(W1)+trace(W2))
    subject to
        intrfPow1 <= epsilon2 - sigma^2;
        intrfPow2 <= epsilon1 - sigma^2;
        desiredPow1 >= epsilon1*(2^rate-1);
        desiredPow2 >= epsilon2*(2^rate-1);
        epsilon1 >= sigma^2;
        epsilon2 >= sigma^2;
cvx_end

if (isnan(W1) | isnan(W2) | isnan(epsilon1) | isnan(epsilon2))
    flag = 1;
else
    flag = 0;
end

end