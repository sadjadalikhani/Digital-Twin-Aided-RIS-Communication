function [phaseCodebook,w2Codebook,w1Codebook] = ...
    codebookGen(angle1,angle2,nBeams,bsAnt,nRIS,nRISbeams)
%% FIND BF VECTOR FOR EACH BS BEAM
deg = linspace(angle1,angle2,nBeams);
bs2risAngle = 50.121;
[~,idx] = min(abs(deg-bs2risAngle));
deg(idx) = bs2risAngle;
f = zeros(bsAnt,nBeams);
cnt = 0;
for theta = deg
    cnt = cnt + 1;
    steerVec = [];
    for n = 0:bsAnt-1
        steerVec = [steerVec; exp(n*1i*2*pi*.5*sinh(deg2rad(theta)))];
    end
    f(:,cnt) = conj(steerVec);
    % PLOT RADIATION PATTERN
    a = zeros(181,1);
    for theta2 = 0:180
        steerVec = [];
        for n = 0:bsAnt-1
            steerVec = [steerVec; exp(n*1i*2*pi*.5*sin(deg2rad(theta2)))];
        end
        a(theta2+1) = abs(steerVec.'*f(:,cnt))^2;
    end
    figure(13)
    polarplot(deg2rad(0:180),a)
    hold on
end
w1Codebook = f;
w2Codebook = f(:,idx);
%% FIND RIS PHASE-SHIFTS FOR THE CENTER OF EACH ZONE 
[~,H2,~,hBSRIS] = RISbeamCSIs(nRISbeams);
psi = zeros(nRIS,nBeams);
for i = 1:nRISbeams
    psi(:,i) = -angle(diag(H2(:,:,i))*hBSRIS*w2Codebook);
end
phaseCodebook = exp(1i.*psi);
%% SAVE CODEBOOKS
% save('phaseCodebook.mat','phaseCodebook')
% save('w2Codebook.mat','w2Codebook')
% save('w1Codebook.mat','w1Codebook')
end
