function [Eprime,flag] = ...
    thetaPrimeOpt(nBS,nRIS,nUE,G1,G2,H3,sigma,rate,epsilon1,epsilon2,W1,W2,H1)
cvx_begin quiet 
    variable Eprime(nRIS+1,nRIS+1) hermitian semidefinite
    expressions intrfPow1 intrfPow2 desiredPow2 intrfPow2beamSweep
    intrfPow1 = 0;
    intrfPow2 = 0;
    desiredPow2 = 0;
    intrfPow2beamSweep = 0;
    for i = 1:nUE
        A1 = double([H1(i,:); diag(G1(i,:))*H3]);
        A2 = double([zeros(1,nBS); diag(G2(i,:))*H3]);
%         A3 = double([H1(i,:); zeros(nRIS,nBS)]);
%         intrfPow1 = intrfPow1 + real(trace(A2*W1*A2'*Eprime));
        intrfPow2 = intrfPow2 + real(trace(A1*W2*A1'*Eprime));
        desiredPow2 = desiredPow2 + real(trace(A2*W2*A2'*Eprime));
%         intrfPow2beamSweep = intrfPow2beamSweep + real(trace(A3*W2*A3'*Eprime));
    end
    subject to
%         intrfPow1 <= epsilon2 - sigma^2;
%         if method == "digitalTwin"
            intrfPow2 <= epsilon1 - sigma^2;
%         elseif method == "beamSweep"
%             intrfPow2beamSweep <= epsilon1 - sigma^2;
%         end
        desiredPow2 >= epsilon2*(2^rate-1);
        for j = 1:nRIS+1
%             Eprime(j,j) <= 1;
            A = zeros(nRIS+1,nRIS+1);
            A(j,j) = 1;
            trace(A*Eprime) == 1;
        end
cvx_end

% Eprime


% Eprime = Eprime./abs(Eprime); 

[Eprime] = rankOneApprox(Eprime,nBS,nRIS,nUE,G1,G2,H3,sigma,rate, ...
                epsilon1,epsilon2,W1,W2,H1,3,3);
Eprime = (Eprime./abs(Eprime));
% eig(Eprime)

if (isnan(Eprime))
    flag = 1;
else
    flag = 0;
end
end