clear all 
close all
tic
%% Packet and transmissions parameters
Header_N_DR8 = 3;       % Header replicas
Code_Rate = 1/3;
Payload = 10;           % Message payload
Header_duration = 0.233; %233 ms long headers
F_duration = 0.05;       %50 ms payload data fragments
Header_ToA_DR8 = Header_N_DR8*Header_duration;
M = 2;             % M = 2 for DR8, M=4 for DR9
Nodes = 1000;
% Simulation_T = 3600; % 1 hour duration
Simulation_T = 30; % 30 segundos
pkct_p_h = 4;      % Packets per hour per end-device
OBW_channels=280;  % No. of OBW channels

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
[Distance, Elevation_Angles, Ground_distance,FootPrint_R]=Satellite_Geometry(H,E);
% X = [1 5 9 13 17 21 25]; %To simulate fewer points
X = [13]; %To simulate fewer points
Distance=Distance(X);

E_angles = [10 20 30 40 50 60 70 80 90];
K_factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
k = sort(interp1(E_angles,K_factor,Elevation_Angles),'descend');

%% Time on air
% ToA_DR8 -> including headers duration
% ToA_DR8_WH -> without headers duration
[ToA_DR8,ToA_DR8_WH] = ToA_Packets_DR8(Payload,Header_ToA_DR8,M); %Função que define o tempo no ar
Transceiver_wait = 6.472/1000; %Tw      
ToA_DR8(1)=ToA_DR8(1) + Transceiver_wait;  % Total on-air time

%% Number of fragments
fragment_duration = 50/1000; %50 ms

fragment_50_ms = floor(ToA_DR8_WH(1)/fragment_duration);  %No. of payload data fragments
% The last fragment may not be equal to 50ms. We need to calculate it.
Last_fragment_duration = ((ToA_DR8_WH(1)/fragment_duration) - fragment_50_ms)*fragment_duration;
fragment_PHY_length  = fragment_50_ms + 1;
fragment_length = Header_N_DR8 + length(Transceiver_wait) + fragment_PHY_length;


%% Simulator
for c=1:1:length(Distance)

    %% Estatisticas
    Success_fragment_received = 0;
    Success_header_received = 0;
    decoded = 0;
    decoded_capture = 0;
    Fragments_received(c) = 0;
    Header_received(c) = 0;

    mu = (1/(Nodes*pkct_p_h)).*Simulation_T;      % inter arrival time
    Inter_arrivals = exprnd(mu,1,Nodes*pkct_p_h); % inter arrival of the traffic in a hour
    First_transmit = rand(1,1)/1000;              % First transmission
    Times = [First_transmit Inter_arrivals];

    %convert inter-arrivals into time stamp 
    TimeStamp = cumsum(Times);                % Time stamp of the traffic in the network
    pack_tx_segments=zeros(length(TimeStamp),fragment_length); % gerando vetor zerado para adicionar os timestamp


    %% Time stamp of the hops (segments) e portadora
    % Definição dos timestamp de cada salto (pack_tx_segments)
    pattern=zeros(length(TimeStamp),fragment_length); % portadoras
    for pack=1:1:length(TimeStamp)   
        for frag = 1:1:fragment_length
                
            % definindo as portadoras para todas as tx
            if frag ~= 1
                dif_track=1;
                while dif_track < 8                       % 8 x 488 = 3.9 kHz spacing
                    pattern(pack,frag) = randi(OBW_channels,1,1);
                    dif_track=abs(pattern(pack,frag)-pattern(pack,frag-1));
                end
            end
    
            % definindo o segmento
            if frag == 1
                pattern(pack,frag) = randi(OBW_channels,1,1);            %First hop of the desired signal
                pack_tx_segments(pack,frag) = TimeStamp (pack);
            elseif frag > 1 && frag <=(Header_N_DR8+1)
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + Header_duration;
            elseif frag == (Header_N_DR8+2)
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + Transceiver_wait;
            else
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + fragment_duration;
            end
        end 
    end
    
    
    %% laço para repetir a comunicação
    for tx=1:1:size(pack_tx_segments,1)  
        %% Vulnerable time
        Transmit = tx;                % Selecting one random device and single transmission instant
        Ts= TimeStamp(Transmit);
            
        % vulnerable time
        Tstart = Ts - ToA_DR8(1);
        Tend = Ts + ToA_DR8(1);
    
        %% Find the number active devices when the desired device was transmitting  
        index_start = find(pack_tx_segments(:,1)>=Tstart & pack_tx_segments(:,1)<=Tend); % Select all the transmission in 2T interval, where T is on-air time
            
        simultaneous = (unique([index_start']));
        target_index = find(simultaneous == Transmit);
        simultaneous(target_index) = [];
        
        %% Collision analysis
        clear seg_simultaneous
        seg_simultaneous = pack_tx_segments(simultaneous,:); % Selecionando intervalos de tempos em que estão sendo transmitidos dados simultaneamente
        target_pattern = pattern(simultaneous,:);            % selecionando as portadoras que estão utilizando os canais
        seg_simultaneous(:,(Header_N_DR8+1))=[];             % removendo TW da simulação
        target_pattern(:,(Header_N_DR8+1))=[];               % removendo TW da simulação
        
        target_collided = zeros(1,size(pack_tx_segments,2));           %Collison counter for Desired signal
        target_discarded = zeros(1,size(pack_tx_segments,2));          %Collison counter for Desired signal for capture effect
     
        % Following Equation (4), (5), (6) to find A_{H}, A_{F} and A_{L}
        for seg=1:1:length(pack_tx_segments(Transmit,:))
        
            clear active_pattern
            clear iscollision
            clear transmission_sim
            clear transmission_sim_frag_device
            clear transmission_sim_header_device
            clear transmission_sim_last_device

            iscollision = 0;
                
            if seg ~= Header_N_DR8+1
                if(seg~=length(pack_tx_segments(Transmit,:)) && seg<=Header_N_DR8)         
                    [row_seg_simultaneous_header, col_seg_simultaneous_header] = find(seg_simultaneous(:,1:Header_N_DR8)>=(pack_tx_segments(Transmit,seg)-Header_duration) ...
                                & seg_simultaneous(:,1:Header_N_DR8)<=(pack_tx_segments(Transmit,seg)+Header_duration));
                                
                    [row_seg_simultaneous_frag, col_seg_simultaneous_frag] = find(seg_simultaneous(:,(Header_N_DR8+1):end-1)>=(pack_tx_segments(Transmit,seg)-F_duration) ...
                                & seg_simultaneous(:,(Header_N_DR8+1):end-1)<=(pack_tx_segments(Transmit,seg)+Header_duration));
            
                    [row_seg_simultaneous_last, col_seg_simultaneous_last] = find(seg_simultaneous(:,end)>=(pack_tx_segments(Transmit,seg)-Last_fragment_duration) ...
                                & seg_simultaneous(:,end)<=(pack_tx_segments(Transmit,seg)+Header_duration));
            
                elseif(seg~=length(pack_tx_segments(Transmit,:)) && seg>Header_N_DR8)
                    [row_seg_simultaneous_header, col_seg_simultaneous_header] = find(seg_simultaneous(:,1:Header_N_DR8)>=(pack_tx_segments(Transmit,seg)-Header_duration) ...
                                & seg_simultaneous(:,1:Header_N_DR8)<=(pack_tx_segments(Transmit,seg)+F_duration));
                                
                    [row_seg_simultaneous_frag, col_seg_simultaneous_frag] = find(seg_simultaneous(:,(Header_N_DR8+1):end-1)>=(pack_tx_segments(Transmit,seg)-F_duration) ...
                                & seg_simultaneous(:,(Header_N_DR8+1):end-1)<=(pack_tx_segments(Transmit,seg)+F_duration));
           
                    [row_seg_simultaneous_last, col_seg_simultaneous_last] = find(seg_simultaneous(:,end)>=(pack_tx_segments(Transmit,seg)-Last_fragment_duration) ...
                                & seg_simultaneous(:,end)<=(pack_tx_segments(Transmit,seg)+F_duration));
              
                else
                    [row_seg_simultaneous_header, col_seg_simultaneous_header] = find(seg_simultaneous(:,1:Header_N_DR8)>=(pack_tx_segments(Transmit,seg)-Header_duration) ...
                                & seg_simultaneous(:,1:Header_N_DR8)<=(pack_tx_segments(Transmit,end))+Last_fragment_duration);
                                
                    [row_seg_simultaneous_frag, col_seg_simultaneous_frag] = find(seg_simultaneous(:,(Header_N_DR8+1):end-1)>=(pack_tx_segments(Transmit,seg)-F_duration) ...
                                & seg_simultaneous(:,(Header_N_DR8+1):end-1)<=(pack_tx_segments(Transmit,end))+Last_fragment_duration);
            
                    [row_seg_simultaneous_last, col_seg_simultaneous_last] = find(seg_simultaneous(:,end)>=pack_tx_segments(Transmit,seg) ...
                                & seg_simultaneous(:,end)<=(pack_tx_segments(Transmit,end))+Last_fragment_duration);
            
                end
                    
                % variável onde estão localizados os fragmentos, este index é referente ao seg_simultaneous
                seg_simultaneous_devices_index = [row_seg_simultaneous_header, col_seg_simultaneous_header; row_seg_simultaneous_frag, col_seg_simultaneous_frag; row_seg_simultaneous_last, col_seg_simultaneous_last];
            
                %% análise de colisões
                if ~isempty(seg_simultaneous_devices_index)
            
                    % nós onde houve a transmissão simultânea
                    simultaneous_device_index = seg_simultaneous_devices_index(:,1);
                    simultaneous_device_unique_index = sort(unique(simultaneous_device_index));
                    simultaneous_device_tx = pack_tx_segments(simultaneous_device_unique_index,:);
            
                    [r_same_pattern, c_same_pattern] = find(target_pattern == pattern(Transmit,seg));
                    same_pattern = [r_same_pattern, c_same_pattern];

                    if(same_pattern)
                        iscollision_index = seg_simultaneous_devices_index(find(ismember(seg_simultaneous_devices_index,same_pattern,'rows')'==1),:);
                        iscollision = sum(ismember(seg_simultaneous_devices_index,same_pattern,'rows')');
                    end
                       
                    if iscollision ~= 0
                        target_collided(seg) = 1;
                        clear Coordinates                    
                        clear rho
                        clear Theta
                                      
                        Coordinates=zeros(iscollision,2);
            
                        rho = sqrt(rand(iscollision,1).*(max(Ground_distance)^2));
                        Theta = rand(iscollision,1)*2*pi;                      
                        Coordinates(:,1) = cos(Theta).*rho;
                        Coordinates(:,2) = sin(Theta).*rho;
            
                        %% Generating location of interfering signals based on uniform distribution of NODES
                        Location_Nodes_Int = sqrt(Coordinates(:,1).^2 + Coordinates(:,2).^2)';
                        dPropogation=zeros(1,iscollision);
                        E_dpro=zeros(1,iscollision);
            
                        % Distance from interfering nodes to satellite
                        for track=1:length(Location_Nodes_Int)
                            dPropogation(1,track) = sqrt(H^2 + Location_Nodes_Int(track).^2);
                            E_dpro(track) = (H*((H+2*R)) - dPropogation(track).^2)./(2.*dPropogation(track).*R);
                        end
            
                        E_AngPro=asind(E_dpro);
                        kC = interp1(E_angles,K_factor,E_AngPro);
            
                        %% Rician fading for interfering signals
                        muC = sqrt(kC./(2.*(kC+1)));    % Mean 
                        sigmaC = sqrt(1./(2.*(kC+1)));  % Variance 
                        hrC=(sigmaC.*randn(1,iscollision)) + muC;
                        hiC=1j.*(sigmaC.*randn(1,iscollision) + muC);
                        h1C=(abs(hrC+hiC)).^2;
            
                        %% total interferance = sum of all the interfering signals 
                        pr_h_g_I = sum(Pt.*h1C.*Gr.*Gt.*((wavelength./(4*pi.*dPropogation)).^eta));
            
                        %% Rician fading for desired signals
                        kD = k(c);    
                        muD = sqrt(kD./(2*(kD+1)));     % Mean 
                        sigmaD = sqrt(1./(2*(kD+1)));   % Variance 
                                
                        hrD=sigmaD*randn(1,1)+muD;
                        hiD=1j.*(sigmaD*randn(1,1)+muD);
                        h1D=(abs(hrD+hiD)).^2;
            
                        %% Received power of desired signal
                        pr_h_g_D = Pt.*h1D.*Gr.*Gt.*((wavelength./(4*pi.*Distance(c))).^eta);
                        if  pr_h_g_D  < (pr_h_g_I*4)
                            %discarded = discarded + 1; %% Destructive collision
                            target_discarded(seg) = 1;
                        end
            
                    end
                end 
            end
        end

        %% Decodificando mensagem
        Success_header = Header_N_DR8 - length(nonzeros(target_collided(1:Header_N_DR8)));       % No. of successfully received headers
        Threshold = size(pack_tx_segments,2) - round(fragment_PHY_length *(1-Code_Rate))-Header_N_DR8 - length(Transceiver_wait);
        Success_fragment = size(target_collided,2) - length(nonzeros(target_collided((Header_N_DR8+2):end)))-Header_N_DR8-1;
        
        % análise somente dos fragmentos
        if(Success_fragment>=Threshold)
            Success_fragment_received = Success_fragment_received+1;
        end
        
        % análise somente dos header
        if (Success_header>=1)
            Success_header_received=Success_header_received+1;
            if(Success_fragment>=Threshold)
                decoded = 1 + decoded;                
            end
        end

        % levando em conta efeito de captura
        Success_header_capture = Header_N_DR8 - length(nonzeros(target_discarded(1:Header_N_DR8)));
        if (Success_header_capture>=1)
            Success_fragment_capture = size(target_discarded,2) - length(nonzeros(target_discarded(Header_N_DR8+2:end)))-Header_N_DR8-1;
            if(Success_fragment_capture>=Threshold)
                decoded_capture = 1 + decoded_capture;
                pack_tx_segments(tx,:) = -100;
            end    
        end

        Fragments_received(c) = Fragments_received(c) + Success_fragment;
        Header_received(c) = Header_received(c) + Success_header;

    end
    %% armazenando estatísticas
    PS_DR8(c)=decoded;   %Simulated overall success probability
    PH_DR8(c)=Success_header_received; %Simulated success probability of headers   
    PF_DR8(c)=Success_fragment_received; %Simulated success probability of data fragments

end

total_fragments = length(TimeStamp)*(fragment_length-3);
total_headers = length(TimeStamp)*(fragment_length-14);

figure(1);
X1 = categorical({'Fragmentos','Header'});
X1 = reordercats(X1,{'Fragmentos','Header'});
Y1 = [Fragments_received, total_fragments; Header_received, total_headers];
graf1 = bar(X1, Y1);
legend({'Recebidos','Enviados'})
title ('Pacotes Independentes - Three Headers (Removendo Transmitidos)');

figure(2)
X2 = categorical({'Fragmentos','Header','Decoded'});
X2 = reordercats(X2,{'Fragmentos','Header','Decoded'});
Y2 = [PF_DR8, total_fragments/14; PH_DR8, total_headers/3; PS_DR8, length(TimeStamp)];
graf2 = bar(X2, Y2);
legend({'Recebidos','Enviados'})
title ('Pacotes Completos - Three Headers (Removendo Transmitidos)');
toc










































