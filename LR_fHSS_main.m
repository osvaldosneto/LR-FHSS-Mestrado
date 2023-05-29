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
F_duration = 0.05;       %50 ms payload data fragments
Header_ToA_DR8 = Header_N_DR8*Header_duration;
Nodes = 150;
Simulation_T = 10; % 1 hour duration
pkct_p_h = 4;      % Packets per hour per end-device
OBW_channels=280;  % No. of OBW channels
M = 2;             % M = 2 for DR8, M=4 for DR9

%% Gains and Pt are converted into linear form

Pt = 10^(14/10)/1000;      % Transmit Power of LoRa 14 dBm
Freq_Band = 868e6;         % 868 MHz (frequency band Europe)
SpeedLight  = 3e8;         % Speed of light
wavelength = SpeedLight/Freq_Band;

Gr=(10.^((22.6)/10));      %22.6 dBi: Gateway at Satellite
Gt=(10.^((2.15)/10));      %2.15 dBi: End-device
eta = 2;

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
ToA_DR8(1)=ToA_DR8(1) + Transceiver_wait;  % Total on-air time

%% Number of fragments
fragment_duration = 50/1000; %50 ms

fragment_50_ms = floor(ToA_DR8_WH(1)/fragment_duration);  %No. of payload data fragments
% The last fragment may not be equal to 50ms. We need to calculate it.
Last_fragment_duration = ((ToA_DR8_WH(1)/fragment_duration) - fragment_50_ms)*fragment_duration;
fragment_PHY_length  = fragment_50_ms + 1;
fragment_length = Header_N_DR8 + length(Transceiver_wait) + fragment_PHY_length;

x = [];
for i=1:1:300
    %% times for transmission
    First_transmit = rand(1,1)/1000;          % First transmission
    mu = (1/(Nodes*pkct_p_h)).*Simulation_T;      % inter arrival time
    Inter_arrivals = exprnd(mu,1,Nodes*pkct_p_h); % inter arrival of the traffic in a hour

    Times = [First_transmit Inter_arrivals];

    %Next step: convert inter-arrivals into time stamp 
    TimeStamp = cumsum(Times);                % Time stamp of the traffic in the network

    [pattern, pr, pack_tx_segments] = Generate_Params(TimeStamp, fragment_length, OBW_channels, ...
    Ground_distance, k, Header_N_DR8, Header_duration, Transceiver_wait, fragment_duration);
    x = [x; Analysys(pattern, pr, pack_tx_segments, Simulation_T, 3, Last_fragment_duration)];
end

disp("final");














