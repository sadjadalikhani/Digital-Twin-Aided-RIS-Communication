%% CONFIGURATION
[h3,h1BS,h1RIS,h2RIS,nBS,nRIS,nUE,sigma,experiments, ...
    rates,bsUsersNum,risUsersNum,locBSusers,locRISusers] = config();
[phaseCodebook,w2Codebook,w1Codebook] = ...
    codebookGen(0,180,16,16,16,16); 
%               angle1,angle2,nBeams,bsAnt,nRIS,nRISbeams
H1_2path = load('./Channels/H1_2path_22Apr2024013001.mat').H1/1e-6;
H1_2path_BSzone = H1_2path(:,:,501*181+1:end);
sumCov = 0;
cumFeas = [];
perfect = 0;
%% MATRICES
[rate1,rate2,txPow1,txPow2,infeasflag,...
    rate1sweep,rate2sweep,txPow1sweep,txPow2sweep,infeasSweep,...
    bsUEidx,risUEidx] = ...
    matrices(length(rates),experiments);  
%% CHANNEL REALIZATIONS
for expr = 1:experiments 
    fprintf('\n\n------------- %g/%g (%g percent) -------------\n\n', ...
        expr,experiments,expr/experiments*100)
    %% CHOOSE TWO RANDOM BS AND RIS USERS
    bsUEidx(expr) = rndRows([1,bsUsersNum],1);
    risUEidx(expr) = rndRows([1,risUsersNum],1);
    %% EXTRACT CHANNELS OF THE CHOSEN TWO USERS
    % 10 paths
    h1 = h1BS(:,:,bsUEidx(expr)); 
    g1 = h1RIS(:,:,bsUEidx(expr)); 
    g2 = h2RIS(:,:,risUEidx(expr)); 
    
    h1_2path = H1_2path_BSzone(:,:,bsUEidx(expr)); 
    sumCov = sumCov + vec(h1_2path-h1)*vec(h1_2path-h1)';
    emprCov = sumCov/expr;
    %% DIGITAL TWIN OPTIMIZATION
    outerflag = 0;
    for rateIdx = 1:length(rates)
        rate = rates(rateIdx);
        %% digital-twin-assisted (10 paths considered)
        w=(unifrnd(-1,1,nBS,1)+1i*unifrnd(-1,1,nBS,1)); 
        W1=w*w'; W2=w*w';
        theta=randn(1,nRIS)+1i*randn(1,nRIS); 
        theta=theta./abs(theta); 
        E=theta'*theta;
        Eprime = [1,theta;theta',E];
        EprimePrev=1e3*ones(nRIS+1,nRIS+1);
        W2Prev=1e1*rand(nBS,nBS)+1i*1e3*rand(nBS,nBS); 
        W1Prev=1e1*rand(nBS,nBS)+1i*1e3*rand(nBS,nBS); 
        cntWhile = 0;
        while (norm(W1-W1Prev,'fro')^2 + norm(W2-W2Prev,'fro')^2 + ...
                norm(Eprime-EprimePrev,'fro')^2 >= 5e-2)
            cntWhile = cntWhile + 1;
            if mod(cntWhile,30) == 0 
                theta=randn(1,nRIS)+1i*randn(1,nRIS); 
                theta=theta./abs(theta); 
                E=theta'*theta;
                Eprime = [1,theta;theta',E];
            elseif cntWhile == 59 
                break
            end
            
            W2Prev = W2;
            W1Prev = W1;
            EprimePrev = Eprime;
            
            if perfect == 1
                [W1,W2,epsilon1,epsilon2,flag] = ...
                    wOpt(nRIS,nBS,nUE,g2,h3,Eprime, ...
                    sigma,h1,rate,g1);
            else
                [W1,W2,epsilon1,epsilon2,flag] = ...
                    wOptRobust(nRIS,nBS,nUE,g2,h3,Eprime, ...
                    sigma,h1,rate,g1,h1_2path,emprCov);
            end
            
            if flag == 1
                W1 = zeros(nBS,nBS);
                W2 = zeros(nBS,nBS);
                outerflag = 1;
                break
            end
            [Eprime,flag] = ...
                thetaPrimeOpt(nBS,nRIS,nUE,g1,g2, ...
                    h3,sigma,rate,epsilon1,epsilon2,W1,W2,h1);
            if flag == 1
                Eprime = [1,zeros(1,nRIS)]'*[1,zeros(1,nRIS)];
                outerflag = 1;
                break
            end
        end
        %% DIGITAL TWIN EVALUATION
        desiredPow1 = trace(h1'*h1*W1);
        desiredPow2 = 0;
        intrf1 = 0;
        intrf2 = 0;
        for i = 1:nUE
            A1 = double([h1(i,:); diag(g1(i,:))*h3]);
            A2 = double([zeros(1,nBS); diag(g2(i,:))*h3]);
            desiredPow2 = desiredPow2 + trace(A2'*Eprime*A2*W2);
            intrf1 = intrf1 + trace(A1'*Eprime*A1*W2);
            intrf2 = intrf2 + trace(A2'*Eprime*A2*W1);
        end
        rate1(rateIdx,expr) = log2(1+real(desiredPow1/(intrf1+sigma^2)));
        rate2(rateIdx,expr) = log2(1+real(desiredPow2/(intrf2+sigma^2)));
        if rate1(rateIdx,expr)<=(rate-.2) | rate2(rateIdx,expr)<=(rate-.2)
            infeasflag(rateIdx,expr) = 1;
        end
        infeasflag(rateIdx,expr) = 1-infeasflag(rateIdx,expr);
        txPow1(rateIdx,expr) = trace(W1);
        if Eprime ~= 0
            txPow2(rateIdx,expr) = trace(W2);
        end
        %% BEAM SWEEPING
        if flag == 0
            [R1sweep,R2sweep,pow1sweep,pow2sweep,infeasflagSweep] = ...
            beamSweep(h1,g1,g2,h3,phaseCodebook,w2Codebook,w1Codebook, ...
                txPow1(rateIdx,expr),txPow2(rateIdx,expr),sigma,nBS,rate);
            rate1sweep(rateIdx,expr) = R1sweep;
            rate2sweep(rateIdx,expr) = R2sweep;
            txPow1sweep(rateIdx,expr) = pow1sweep;
            txPow2sweep(rateIdx,expr) = pow2sweep;
            infeasSweep(rateIdx,expr) = 1-infeasflagSweep;
        end
    end
    %% DIGITAL TWIN RESULTS FOR THIS REALIZATION
    fprintf("\n\n\n DT PERFORMANCE- CURRENT REALIZATION:\n")  
    fprintf("\n user 1 rate (bit/s/Hz): "); 
    fprintf("%g, ",rate1(:,expr)); 
    fprintf("\n user 2 rate (bit/s/Hz): "); 
    fprintf("%g, ",rate2(:,expr)); 
    fprintf("\n user 1 tx power (dBm): "); 
    fprintf("%g, ",10*log10(1e3*txPow1(:,expr))); 
    fprintf("\n user 2 tx power (dBm): "); 
    fprintf("%g, ",10*log10(1e3*txPow2(:,expr)));
    fprintf("\n total tx power (dBm): "); 
    fprintf("%g, ",10*log10(1e3*(txPow1(:,expr)+txPow2(:,expr))));
    fprintf("\n feasibility rate: "); 
    fprintf("%g, ",infeasflag(:,expr)); 
    fprintf("\n cumulative feasibility rate: "); 
    fprintf("%g, ",nanmean(infeasflag(:,1:expr),2)); 

    %% BEAM SWEEPING RESULTS FOR THIS REALIZATION
    fprintf("\n\n\n BEAM SWEEPING PERFORMANCE- CURRENT REALIZATION:\n")
    fprintf("\n user 1 rate (bit/s/Hz): "); 
    fprintf("%g, ",rate1sweep(:,expr)); 
    fprintf("\n user 2 rate (bit/s/Hz): "); 
    fprintf("%g, ",rate2sweep(:,expr)); 
    fprintf("\n user 1 tx power (dBm): "); 
    fprintf("%g, ",10*log10(1e3*txPow1sweep(:,expr))); 
    fprintf("\n user 2 tx power (dBm): "); 
    fprintf("%g, ",10*log10(1e3*txPow2sweep(:,expr)));
    fprintf("\n total tx power (dBm): "); 
    fprintf("%g, ",10*log10(1e3*(txPow1sweep(:,expr)+txPow2sweep(:,expr))));
    fprintf("\n feasibility rate: "); 
    fprintf("%g, ",infeasSweep(:,expr)); 

end

save('rate1.mat','rate1')
save('rate2.mat','rate2')
save('txPow1.mat','txPow1')
save('txPow2.mat','txPow2')
save('infeasflag.mat','infeasflag')

save('rate1sweep.mat','rate1sweep')
save('rate2sweep.mat','rate2sweep')
save('txPow1sweep.mat','txPow1sweep')
save('txPow2sweep.mat','txPow2sweep')
save('infeasSweep.mat','infeasSweep')

%% Bar Plot

% First load the perfect design rates and then robust and beam sweeping 
if perfect == 1
rate1_opt = rate1; 
rate2_opt = rate2;
else
rate1_rob = rate1;
rate2_rob = rate2;
rate1_bs = rate1sweep;
rate2_bs = rate2sweep;
end

% Then run this part 
x = [1,1.5,2,2.5,3,3.5,4];
y = [nanmean(min(rate1_opt,x.')+min(rate2_opt,x.'),2), ...
     nanmean(min(rate1_rob,x.')+min(rate2_rob,x.'),2), ...
     nanmean(min(rate1_bs,x.')+min(rate2_bs,x.'),2)];

barPlot(x, y)

%% CDF Plot

% First load the perfect design rates and then robust and beam sweeping 
if perfect == 1
rate1_opt = rate1(1,:); 
rate2_opt = rate2(1,:);
else
rate1_rob = rate1(1,:);
rate2_rob = rate2(1,:);
rate1_bs = rate1sweep(1,:);
rate2_bs = rate2sweep(1,:);
end

% Then run this part 
X1 = cdfplot(min(rate1_opt,1)+min(rate2_opt,1)).XData;
Y1 = cdfplot(min(rate1_opt,1)+min(rate2_opt,1)).YData;

X2 = cdfplot(min(rate1_rob,1)+min(rate2_rob,1)).XData;
Y2 = cdfplot(min(rate1_rob,1)+min(rate2_rob,1)).YData;

X3 = cdfplot(min(rate1_bs,1)+min(rate2_bs,1)).XData;
Y3 = cdfplot(min(rate1_bs,1)+min(rate2_bs,1)).YData;

cdf(X1, Y1, X2, Y2, X3, Y3)


