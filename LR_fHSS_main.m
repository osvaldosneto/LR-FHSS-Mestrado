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
Nodes = 1000;
Simulation_T = 3600;   % 1 hour duration
pkct_p_h = 4;      % Packets per hour per end-device
OBW_channels=280;  % No. of OBW channels
M = 2;             % M = 2 for DR8, M=4 for DR9
t_hour = 30;

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


n_sim = 10;
sim = zeros(n_sim,Nodes*pkct_p_h);
for i=1:n_sim

    %% times for transmission
%     mu = 3600/(pkct_p_h);
    mu = 30/(pkct_p_h);
%     Inter_arrivals = exprnd(mu,Nodes,round((pkct_p_h*Simulation_T)/3600)); % inter arrival of the traffic in a hour
    Inter_arrivals = exprnd(mu,Nodes,pkct_p_h); % inter arrival of the traffic in a hour
    Times = Inter_arrivals;
    Times(:,2:end) = Times(:,2:end)+time_tx_pack;

    %Next step: convert inter-arrivals into time stamp 
    TimeStamp = cumsum(Times,2);                % Time stamp of the traffic in the network

    [pattern, pr, pack_tx_segments, distance] = Generate_Params(TimeStamp, fragment_length, OBW_channels, ...
    Ground_distance, k, Header_N_DR8, Header_duration, Transceiver_wait, fragment_duration, ...
    Nodes,pkct_p_h);
    a = Analysys(pattern, pr, pack_tx_segments, Simulation_T, 3, Last_fragment_duration);
    sim(i,1:length(a)) = a';
end
toc

%% extraindo estatísticas da transmissão
d_slice = [1.1e6 1.3e6 1.5e6 1.7e6 1.9e6 2.1e6 2.3e6];
dif1 = 0; devices1 = 0;
dif2 = 0; devices2 = 0;
dif3 = 0; devices3 = 0;
dif4 = 0; devices4 = 0;
dif5 = 0; devices5 = 0;
dif6 = 0; devices6 = 0;
dif7 = 0; devices7 = 0;
dif8 = 0; devices8 = 0;

for n=1:1:size(sim,1)
    for d=1:1:length(d_slice)

        switch d
            case 1
                devices = find(distance < d_slice(d));
                dif1 = dif1 + length(devices)-length(setdiff(devices,sim(n,:)));
                devices1 = devices1 + length(devices);

            case 2
                devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                dif2 = dif2 + length(devices)-length(setdiff(devices,sim(n,:)));
                devices2 = devices2 + length(devices);

            case 3
                devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                dif3 = dif3 + length(devices)-length(setdiff(devices,sim(n,:)));
                devices3 = devices3 + length(devices);

            case 4
                devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                dif4 = dif4 + length(devices)-length(setdiff(devices,sim(n,:)));
                devices4 = devices4 + length(devices);

            case 5
                devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                dif5 = dif5 + length(devices)-length(setdiff(devices,sim(n,:)));
                devices5 = devices5 + length(devices);

            case 6
                devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                dif6 = dif6 + length(devices)-length(setdiff(devices,sim(n,:)));
                devices6 = devices6 + length(devices);

            case 7
                devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                dif7 = dif7 + length(devices)-length(setdiff(devices,sim(n,:)));
                devices7 = devices7 + length(devices);

            case 8
                devices = find(distance >= d_slice(d));
                dif8 = dif8 + length(devices)-length(setdiff(devices,sim(n,:)));
                devices8 = devices8 + length(devices);
        end

    end

end


dif1 = dif1/devices1;
dif2 = dif2/devices2;
dif3 = dif3/devices3;
dif4 = dif4/devices4;
dif5 = dif5/devices5;
dif6 = dif6/devices6;
dif7 = dif7/devices7;
dif8 = dif8/devices8;

d_plot = [1e6 1.2e6 1.4e6 1.6e6 1.8e6 2.0e6];
prob = [dif1 dif2 dif3 dif5 dif6 dif7];

figure(1)
plot(d_plot, prob, 'bo')
axis([1e6 2.0e6 0 1])
grid on
ylabel('Success probability', 'Interpreter', 'Latex');
xlabel('Distance from node to satellite (km)', 'Interpreter', 'Latex');



