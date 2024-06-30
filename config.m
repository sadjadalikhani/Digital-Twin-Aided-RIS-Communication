function [hBSRIS,h1BS,h1RIS,h2RIS,nBS,nRIS,nUE,sigma,experiments, ...
    rates,bsUsersNum,risUsersNum,locBSusers,locRISusers] = config()

% % h[ZONE][SOURCE][#PATHS] %%%%%%%%% ZONE(BS USERS:1,RIS USERS:2)
% h1BS1 = (load('H1_1pathBSusers_04Apr2024181839.mat').H1)./1e-6;
% h1RIS1 = (load('H2_1pathBSusers_04Apr2024181839.mat').H2)./1e-6;
% h2RIS1 = (load('H2_1pathRISusers_04Apr2024230745.mat').H2)./1e-6; 
% h1BS10 = (load('H1_10pathsBSusers_04Apr2024185827.mat').H1)./1e-6;
% h1RIS10 = (load('H2_10pathsBSusers_04Apr2024185827.mat').H2)./1e-6;
% h2RIS10 = (load('H2_10pathsRISusers_04Apr2024231047.mat').H2)./1e-6;
% hBSRIS1 = (load('hBSBS_1path_05Apr2024003050.mat').hBSBS)./1e-6;
% hBSRIS10 = (load('hBSBS_10paths_05Apr2024003143.mat').hBSBS)./1e-6;
% locBSusers = load('loc_1pathBSusers_04Apr2024181839.mat').loc;
% locRISusers = load('loc_1pathRISusers_04Apr2024230745.mat').loc;


H1 = (load('./Channels/H1_10path_17Apr2024024034.mat').H1)./1e-6;
H2 = (load('./Channels/H2_10path_17Apr2024024034.mat').H2)./1e-6;
hBSRIS = (load('./Channels/hBSBS_10path_17Apr2024024034.mat').hBSBS)./1e-6;
locs = load('./Channels/loc_10path_17Apr2024024034.mat').loc;

h1BS = H1(:,:,501*181+1:end);
h1RIS = H1(:,:,1:501*181);
h2RIS = H2(:,:,1:501*181);
locBSusers = locs(501*181+1:end,:);
locRISusers = locs(1:501*181,:);

bsUsersNum = size(h1RIS,3);
risUsersNum = size(h2RIS,3);


nBS = 16;
nRIS = 16;
nUE = 1;

sigma = 1e-2;
experiments = 2000;
% rates = [1,1.5,2,2.5,3,3.5,4];
rates = [1];

end