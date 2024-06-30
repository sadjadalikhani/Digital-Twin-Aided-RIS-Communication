function [b] = gaussRand(a)

    [eigVec,eigVal] = eig(a);
    r = (randn(size(a,1),1) + 1i*randn(size(a,1),1))/sqrt(2);
    b = eigVec*eigVal^(1/2)*(r*r')*eigVal^(1/2)'*eigVec';

end