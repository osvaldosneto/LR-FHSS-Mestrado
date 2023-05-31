clear all 
close all

clear TimeStamp
clear pack_tx_segments

%% Packet and transmissions parameters
tic
Payload = 10;           % Message payload
Header_N_DR8 = 3;       % Header replicas
Code_Rate = 1/3;
Header_duration = 0.233; %233 ms long headers
fragment_duration = 50/1000; %50 ms
Header_ToA_DR8 = Header_N_DR8*Header_duration;
Nodes = 150;
Simulation_T = 3600; % 1 hour duration
pkct_p_h = 4;      % Packets per hour per end-device
OBW_channels=280;  % No. of OBW channels
M = 2;             % M = 2 for DR8, M=4 for DR9

%% parameters for satellites
E = 10:1:90;               %Elevation Angles
R = 6378e3;                % Radius of earth
H = 780e3;                 %Orbital height  

%% Distance from user to satellite as function of elevation angle
[Distance, Elevation_Angles, Ground_distance, FootPrint_R]=Satellite_Geometry(H,E);

E_angles = [10 20 30 40 50 60 70 80 90];
K_factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
k = sort(interp1(E_angles,K_factor,Elevation_Angles),'descend');

%% Time on air
[ToA_DR8,ToA_DR8_WH] = ToA_Packets_DR8(Payload,Header_ToA_DR8,M); 
% ToA_DR8 -> including headers duration
% ToA_DR8_WH -> without headers duration
Transceiver_wait = 6.472/1000; %Tw      

%% Number of fragments
fragment_50_ms = floor(ToA_DR8_WH(1)/fragment_duration);  %No. of payload data fragments
% The last fragment may not be equal to 50ms. We need to calculate it.
Last_fragment_duration = ((ToA_DR8_WH(1)/fragment_duration) - fragment_50_ms)*fragment_duration;
fragment_PHY_length  = fragment_50_ms + 1;
fragment_length = Header_N_DR8 + length(Transceiver_wait) + fragment_PHY_length;
time_tx_pack = Header_N_DR8*Header_duration + Transceiver_wait + fragment_PHY_length*fragment_duration;


n_sim = 50;
sim = ones(n_sim,1).*NaN;
% parfor i=1:n_sim

    %% times for transmission
%     mu = (1/(Nodes*pkct_p_h)).*Simulation_T;      % inter arrival time
    mu = rand();
    Inter_arrivals = exprnd(mu,Nodes,pkct_p_h); % inter arrival of the traffic in a hour
    Times = Inter_arrivals;
    Times(:,2:end) = Times(:,2:end)+time_tx_pack;

    %Next step: convert inter-arrivals into time stamp 
    TimeStamp = cumsum(Times,2);                % Time stamp of the traffic in the network

    [pattern, pr, pack_tx_segments, distance] = Generate_Params(TimeStamp, fragment_length, OBW_channels, ...
    Ground_distance, k, Header_N_DR8, Header_duration, Transceiver_wait, fragment_duration, ...
    Nodes,pkct_p_h);
    a = Analysys(pattern, pr, pack_tx_segments, Simulation_T, 3, Last_fragment_duration);
%     sim(i) = length(a);
% end
toc
