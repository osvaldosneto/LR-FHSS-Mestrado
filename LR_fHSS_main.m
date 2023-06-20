clear all 
close all

clear TimeStamp
clear pack_tx_segments

%% Packet and transmissions parameters
tic
Payload = 10;           % Message payload
Header_N_DR8 = 3;       % Header replicas
Nodes = 2500;           % número de nós na rede
pkct_p_h = 4;           % Packets per hour per end-device
OBW_channels=280;       % No. of OBW channels
Simulation_T = 60;      % duração em segundos

n_sim = 10;
sim = zeros(n_sim,Nodes*pkct_p_h);
sim_SIR = zeros(n_sim,Nodes*pkct_p_h);

%% simulando
for i=1:n_sim
    [pattern, pr, pack_tx_segments, distance, last_fragment_duration] = Generate_Params( ...
        OBW_channels, Header_N_DR8, Nodes, pkct_p_h, Payload, Simulation_T);

    analysys = Analysys(pattern, pr, pack_tx_segments, Header_N_DR8, last_fragment_duration, 0);
    sim(i,1:length(analysys)) = analysys';

    analysys_SIR = Analysys(pattern, pr, pack_tx_segments, Header_N_DR8, last_fragment_duration, 1);
    sim_SIR(i,1:length(analysys_SIR)) = analysys_SIR';
end

[d_plot, prob] = Generate_Statistics(sim, distance);
[d_plot_SIR, prob_SIR] = Generate_Statistics(sim_SIR, distance);

figure(1);
plot(d_plot, prob, 'bo');
axis([1e6 2.0e6 0 1]);
grid on;
hold on;
plot(d_plot_SIR, prob_SIR, 'ro');
legend('Sem SIR', 'Com SIR');
ylabel('Success probability', 'Interpreter', 'Latex');
xlabel('Distance from node to satellite (km)', 'Interpreter', 'Latex');

toc