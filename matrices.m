function [rate1,rate2,txPow1,txPow2,infeasflag,...
    rate1sweep,rate2sweep,txPow1sweep,txPow2sweep,infeasSweep,...
    bsUEidx,risUEidx] = matrices(ratesNum,experiments)

rate1 = zeros(ratesNum,experiments);
rate2 = zeros(ratesNum,experiments);
txPow1 = zeros(ratesNum,experiments);
txPow2 = zeros(ratesNum,experiments);
infeasflag = zeros(ratesNum,experiments);

rate1sweep = NaN*zeros(ratesNum,experiments);
rate2sweep = NaN*zeros(ratesNum,experiments);
txPow1sweep = NaN*zeros(ratesNum,experiments);
txPow2sweep = NaN*zeros(ratesNum,experiments);
infeasSweep = NaN*zeros(ratesNum,experiments);

bsUEidx = zeros(experiments,1);
risUEidx = zeros(experiments,1);
end