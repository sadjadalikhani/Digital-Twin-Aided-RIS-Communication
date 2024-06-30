function [R1_opt,R2_opt,txPow1_opt,txPow2_opt,infeasflag] = ...
    beamSweep(h1,g1,g2,h3,phaseCodebook,w2Codebook,w1Codebook, ...
        txPow1DT,txPow2DT,sigma,nBS,targetRate)
w1CodebookSize = size(w1Codebook,2);
phaseCodebookSize = size(phaseCodebook,2);
%% MATRIX DEFINITIONS

R1 = zeros(phaseCodebookSize,w1CodebookSize);
R2 = zeros(phaseCodebookSize,w1CodebookSize);
txPow1 = zeros(phaseCodebookSize,w1CodebookSize);
txPow2 = zeros(phaseCodebookSize,w1CodebookSize);

%% EXHAUSTIVE SEARCH OVER BF VECTORS AND PHASE-SHIFTS

for i = 1:w1CodebookSize
    for j = 1:phaseCodebookSize

        w1 = sqrt((txPow1DT+txPow2DT)/2)*w1Codebook(:,i)/norm(w1Codebook(:,i));
        w2 = sqrt((txPow1DT+txPow2DT)/2)*w2Codebook/norm(w2Codebook);

        theta = phaseCodebook(:,j);

        Eprime = [1,theta.';conj(theta),theta*theta'];
        W1 = w1*w1';
        W2 = w2*w2';

        desiredPow1 = trace(h1'*h1*W1);
        A1 = double([h1; diag(g1)*h3]);
        A2 = double([zeros(1,nBS); diag(g2)*h3]);
        desiredPow2 = trace(A2'*Eprime*A2*W2);
        intrf1 = trace(A1'*Eprime*A1*W2);
        intrf2 = trace(A2'*Eprime*A2*W1);
        
        R1(j,i) = log2(1+real(desiredPow1/(intrf1+sigma^2)));
        R2(j,i) = log2(1+real(desiredPow2/(intrf2+sigma^2)));

        txPow1(j,i) = trace(W1); 
        txPow2(j,i) = trace(W2); 
    end
end

%% FIND THE SATISFIED RATE WITH THE MINIMUM POWER

R1_1 = sign((R1-targetRate));
R1_1(R1_1<0)=0;
R2_1 = sign((R2-targetRate));
R2_1(R2_1<0)=0;
[idx1,idx2] = find(R1_1.*R2_1==1);
infeasflag = 0;
if isempty([idx1,idx2])
    infeasflag = 1;
    sumSquaredR = R1+R2;
    [idx3,idx4] = find(sumSquaredR==max(sumSquaredR(:)));
    R1_opt = R1(idx3(1),idx4(1));
    R2_opt = R2(idx3(1),idx4(1));
    txPow1_opt = txPow1(idx3(1),idx4(1));
    txPow2_opt = txPow2(idx3(1),idx4(1));
else
    txPowMin = zeros(length(idx1),1);
    for q = 1:length(idx1)
        txPowMin(q) = txPow1(idx1(q),idx2(q))+txPow2(idx1(q),idx2(q));
    end
    [~,txPowMinIdx] = min(txPowMin);
    idx1(txPowMinIdx)
    idx2(txPowMinIdx)
    R1_opt = R1(idx1(txPowMinIdx),idx2(txPowMinIdx));
    R2_opt = R2(idx1(txPowMinIdx),idx2(txPowMinIdx));
    txPow1_opt = txPow1(idx1(txPowMinIdx),idx2(txPowMinIdx));
    txPow2_opt = txPow2(idx1(txPowMinIdx),idx2(txPowMinIdx));
end

%%
end
