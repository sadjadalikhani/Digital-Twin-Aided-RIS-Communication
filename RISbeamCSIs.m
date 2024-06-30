function [H1,H2,loc,hBSRIS] = RISbeamCSIs(risZonePartitions)
%% PARTITION RIS ZONE AND FIND THE USER AT THE CENTER OF EACH ZONE
% VERTICAL
a = linspace(1,181,4+1);
b = [];
for i = 1:4
b = [b floor((a(i+1)-a(i))/2+a(i))];
end
% HORIZONTAL
c = linspace(900,1400,risZonePartitions/4+1);
d = [];
for i = 1:risZonePartitions/4
d = [d floor((c(i+1)-c(i))/2+c(i))];
end
%% COMPUTE CORRESPONDING CHANNELS
addpath('DeepMIMO_functions')
format shortG
r = d;
num = 16;
dSize = 1;
APSize = 16;
H1 = zeros(dSize,APSize,num);
H2 = zeros(dSize,APSize,num);
loc = zeros(num,2);

cnt3 = 0; 
cnt = -1;
for i = r
    cnt = cnt + 1;
    (cnt+1)/length(r)*100
    dataset_params = read_params('parameters.m');
    dataset_params.active_user_first = i;
    dataset_params.active_user_last = i;
    
    [DeepMIMO_dataset, ~] = DeepMIMO_generator(dataset_params);

    for j = b
        cnt3 = cnt3 + 1;
        H1(:,:,cnt3) = DeepMIMO_dataset{1}.user{j}.channel;
        H2(:,:,cnt3) = DeepMIMO_dataset{2}.user{j}.channel;
        pos = DeepMIMO_dataset{1}.user{j}.loc;
        loc(cnt3,:) = pos(1:2);
    end
end
hBSRIS = DeepMIMO_dataset{1}.basestation{2}.channel;
end
%%
% -------------------------- Output Examples -----------------------------%
% DeepMIMO_dataset{i}.user{j}.channel % Channel between BS i - User j
% %  (# of User antennas) x (# of BS antennas) x (# of OFDM subcarriers)
%
% DeepMIMO_dataset{i}.user{j}.params % Parameters of the channel (paths)
% DeepMIMO_dataset{i}.user{j}.LoS_status % Indicator of LoS path existence
% %     | 1: LoS exists | 0: NLoS only | -1: No paths (Blockage)|
%
% DeepMIMO_dataset{i}.user{j}.loc % Location of User j
% DeepMIMO_dataset{i}.loc % Location of BS i
%
% % BS-BS channels are generated only if (params.enable_BSchannels == 1):
% DeepMIMO_dataset{i}.basestation{j}.channel % Channel between BS i - BS j
% DeepMIMO_dataset{i}.basestation{j}.loc
% DeepMIMO_dataset{i}.basestation{j}.LoS_status
%
% % Recall that the size of the channel vector was given by 
% % (# of User antennas) x (# of BS antennas) x (# of OFDM subcarriers)
% % Each of the first two channel matrix dimensions follows a certain 
% % reshaping sequence that can be obtained by the following
% % 'antennamap' vector: Each entry is 3 integers in the form of 
% % 'xyz' where each representing the antenna number in x, y, z directions
% antennamap = antenna_channel_map(params.num_ant_BS(1), ...
%                                  params.num_ant_BS(2), ...
%                                  params.num_ant_BS(3), 1);
%
% -------------------------- Dynamic Scenario ----------------------------%
%
% DeepMIMO_dataset{f}{i}.user{j}.channel % Scene f - BS i - User j
% % Every other command applies as before with the addition of scene ID
% params{f} % Parameters of Scene f
%
% ------------------------------------------------------------------------%