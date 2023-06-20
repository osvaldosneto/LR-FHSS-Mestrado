clear all 
close all

clear TimeStamp
clear pack_tx_segments

%% Packet and transmissions parameters
tic
Payload = 10;           % Message payload
Header_N_DR8 = 3;       % Header replicas
Nodes = (200:100:1600); % número de nós na rede
n_sim = 5;              % número de simulações
pkct_p_h = 4;           % Packets per hour per end-device
OBW_channels=280;       % No. of OBW channels
Simulation_T = 60;      % duração em segundos

prob = zeros(1, length(Nodes));
prob_SIR = zeros(1, length(Nodes));
n_colisions_array = zeros(1, length(Nodes));
total_devices_tx = zeros(1, length(Nodes));

for n=1:length(Nodes)

    received_pakage = 0;
    received_pakage_SIR = 0;

    n_colisions = 0;
    n_colisions_SIR = 0;
    
    %% simulando
    for i=1:n_sim
        [pattern, pr, pack_tx_segments, distance, last_fragment_duration] = Generate_Params( ...
            OBW_channels, Header_N_DR8, Nodes(n), pkct_p_h, Payload, Simulation_T);
    
        [received, colisions] = Analysys(pattern, pr, pack_tx_segments, Header_N_DR8, last_fragment_duration, 0);
        received_pakage = received_pakage + length(received);
        n_colisions = n_colisions + colisions;
    
        [received_SIR, colisions_SIR] = Analysys(pattern, pr, pack_tx_segments, Header_N_DR8, last_fragment_duration, 1);
        received_pakage_SIR = received_pakage_SIR + length(received_SIR);
    end

    prob(n) = received_pakage/(n_sim*length(pack_tx_segments));
    prob_SIR(n) = received_pakage_SIR/(n_sim*length(pack_tx_segments));
    n_colisions_array(n) = n_colisions;
    total_devices_tx(n) = n_sim*length(pack_tx_segments);

end

%% plot de estatísticas
figure(1);
plot(Nodes, prob, 'bo');
axis([(Nodes(1)-50) Nodes(length(Nodes))+50 0 1]);
grid on;
hold on;
plot(Nodes, prob_SIR, 'ro');
legend('Normal', 'SIR');
ylabel('Success probability', 'Interpreter', 'Latex');
xlabel('Number of devices', 'Interpreter', 'Latex');
title('SIMULATION LRFHSS');

figure(2);
bar(Nodes, n_colisions_array);
ylabel('Number of Collisions', 'Interpreter', 'Latex');
xlabel('Number of devices', 'Interpreter', 'Latex');
title('HISTOGRAM - COLLISIONS LRFHSS');

toc
