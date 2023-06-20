% Função para gerar as portadoras, distâncias, potência recebida e o
% coeficiente de desvanecimento para cada transmissão de fragmento

function [pattern, pr_h_g_D, pack_tx_segments, distance, last_fragment_duration, time_tx_pack] = Generate_Params( ...
    OBW_channels, header_size, nodes, pkct_p_h, payload, simulation_time)

    % variáveis fixas
    Gr=(10.^((22.6)/10));           % 22.6 dBi: Gateway at Satellite
    Gt=(10.^((2.15)/10));           % 2.15 dBi: End-device
    R = 6378e3;                     % Radius of earth
    H = 780e3;                      % Orbital height
    E = 10:1:90;                    % Elevation Angles
    Pt = 10^(14/10)/1000;           % Transmit Power of LoRa 14 dBm
    Freq_Band = 868e6;              % 868 MHz (frequency band Europe)
    SpeedLight  = 3e8;              % Speed of light
    eta = 2;                        % expoente de atenuação do caminho
    fragment_duration = 50/1000;    % 50 ms
    Header_duration = 0.233;        % 233 ms long headers
    M = 2;                          % M = 2 for DR8, M=4 for DR9

    wavelength = SpeedLight/Freq_Band;
    Header_ToA_DR8 = header_size*Header_duration;

    %% Distance from user to satellite as function of elevation angle
    [Distance, Elevation_Angles, Ground_distance, FootPrint_R]=Satellite_Geometry(H,E);
    E_angles = [10 20 30 40 50 60 70 80 90];
    K_factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
    k = sort(interp1(E_angles,K_factor,Elevation_Angles),'descend');

    %% Time on air
    [ToA,ToA_WH] = ToA_Packets_DR8(payload,Header_ToA_DR8,M); 
    Transceiver_wait = 6.472/1000; %Tw

    fragment_50_ms = floor(ToA_WH(1)/fragment_duration);  %No. of payload data fragments
    fragment_PHY_length  = fragment_50_ms + 1;
    last_fragment_duration = ((ToA_WH(1)/fragment_duration) - fragment_50_ms)*fragment_duration;
    fragment_length = header_size + length(Transceiver_wait) + fragment_PHY_length;

    time_tx_pack = header_size*Header_duration + Transceiver_wait + fragment_PHY_length*fragment_duration;

    %% definimos as variáveis a utilizar na simulação
    pattern = zeros(nodes*pkct_p_h,fragment_length);
    pr_h_g_D = zeros(nodes*pkct_p_h,fragment_length);
    pack_tx_segments=zeros(nodes*pkct_p_h,fragment_length);
    distance = zeros(nodes*pkct_p_h,1);

    %% cálculo distância entre os dispositivos e satélite
    for d=1:pkct_p_h:(nodes*pkct_p_h-1)
        rho = sqrt(rand(1,1)*(max(Ground_distance)^2));
        Theta = rand(1,1)*2*pi;                      
        Coordinates(:,1) = cos(Theta).*rho;
        Coordinates(:,2) = sin(Theta).*rho;
        distance(d:d+pkct_p_h-1,1) = sqrt(Coordinates(:,1)^2 + Coordinates(:,2)^2)';
    end

    %% definindo timestamp do inicio de cada transmissão

    % utilizando distribuição randômica para gerar pacotes
       time_stamp = (simulation_time-time_tx_pack)*rand(pkct_p_h,nodes);
       time_stamp = sort(time_stamp);
       
       [l, c] = find(abs(diff(time_stamp)) < time_tx_pack);
       tx_same_time = unique(c);
   
       for i=1:1:length(tx_same_time)
           
           while ~isempty(c)
               time_stamp(:,tx_same_time(i)) = (simulation_time-time_tx_pack)*rand(pkct_p_h,1);
               time_stamp(:,tx_same_time(i)) = sort(time_stamp(:,tx_same_time(i)));
               [l, c] = find(abs(diff(time_stamp(:,tx_same_time(i)))) < time_tx_pack);
              
           end
           [l, c] = find(abs(diff(time_stamp)) < time_tx_pack);
       end

    % utilizando distribuição exponencial para gerar pacotes
%       mu = 30/(pkct_p_h);
%       Inter_arrivals = exprnd(mu,nodes,pkct_p_h); % inter arrival of the traffic in a hour
%       Times = Inter_arrivals;
%       Times(:,2:end) = Times(:,2:end)+time_tx_pack;
%       time_stamp = cumsum(Times,2);
    
    for pack=1:1:length(pack_tx_segments)
        
        %% definindo timstamp para cada transmissão
        for frag = 1:1:fragment_length
            if frag == 1
                pack_tx_segments(pack,frag) = time_stamp(pack);
            elseif frag > 1 && frag <=(header_size+1)
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + Header_duration;
                   
            elseif frag == (header_size+2)
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + Transceiver_wait;
    
            else
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + fragment_duration;
            end
    
            %% definição das portadoras em cada tx
            if frag == 1
                pattern(pack,frag) = randi(OBW_channels,1,1);
            elseif frag == (header_size+1)
                pattern(pack,frag) = 0; % não exist portadora twait
            else
                dif_track=0;
                 while dif_track < 8  % 8 x 488 = 3.9 kHz spacing
                    pattern(pack,frag) = randi(OBW_channels,1,1);
                     dif_track=abs(pattern(pack,frag)-pattern(pack,frag-1));
                 end
            end
    
            % potência recebida para cada segmento
            kD = k(randi(length(k))); % fator Riciano
            muD = sqrt(kD/(2*(kD+1)));          % Mean 
            sigmaD = sqrt(1/(2*(kD+1)));        % Variance                          
            hrD=sigmaD*randn(1,1)+muD;          % parte real coeficiente de desvanecimento
            hiD=1j*(sigmaD*randn(1,1)+muD);     % aprte imaginária coeficiente de desvanecimento
            h1D=(abs(hrD+hiD)).^2;    % módulo do coeficiente de desvanecimento Riciano
            pr_h_g_D(pack,frag) = Pt*h1D*Gr*Gt*((wavelength/(4*pi*distance(pack)))^eta);

        end

    end

end
