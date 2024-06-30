function [F] = rankOneApprox(F,nBS,nRIS,nUE,G1,G2,H3,sigma,rate, ...
                epsilon1,epsilon2,W1,W2,H1,constrNum,label)
            
    [a, b] = size(F);
    if ~isnan(F)
        eigF = eig(F);
        flag2 = 0;
        cnt = 0;
        flag1 = 0;
        if length(eigF(abs(eigF)<1e-6 == 0)) > 1
            [eigVec,eigVal] = eig(F);
            [maxVal,idx] = max(diag(eigVal));
            rankOneF = maxVal*eigVec(:,idx)*eigVec(:,idx)';
%             xxx = ConstrCheck(F,nBS,nRIS,nUE,G1,G2,H3,sigma,rate, ...
%                 epsilon1,epsilon2,W1,W2,H1,label);
%             if xxx == ones(constrNum,1)
%                 flag1 = 1;
%             end
% 
%             if flag1 == 0
%                 while flag2 == 0
%                     cnt = cnt + 1;
%                     if cnt == 20
%                         rankOneF = zeros(a,b);
%                         break
%                     end
%                     [rankOneF] = gaussRand(F);
%                     xxx = ConstrCheck(rankOneF,nBS,nRIS,nUE,G1,G2,H3,sigma,rate, ...
%                         epsilon1,epsilon2,W1,W2,H1,label);
%                     if xxx == ones(constrNum,1) 
%                         flag2 = 1;
%                     end
%                 end
%             end
            F = rankOneF;
        end
    else
        F = zeros(a,b);
    end
end