% Função para gerar as portadoras, distâncias, potência recebida e o
% coeficiente de desvanecimento para cada transmissão de fragmento

function [pattern, pr_h_g_D, pack_tx_segments] = Generate_Params(TimeStamp, ...
    fragment_length, OBW_channels, Ground_distance, k, header_size, Header_duration, ...
    Transceiver_wait, fragment_duration)

    % variáveis fixas
    Gr=(10.^((22.6)/10));      % 22.6 dBi: Gateway at Satellite
    Gt=(10.^((2.15)/10));      % 2.15 dBi: End-device
    R = 6378e3;                % Radius of earth
    H = 780e3;                 % Orbital height  
    Pt = 10^(14/10)/1000;      % Transmit Power of LoRa 14 dBm
    Freq_Band = 868e6;         % 868 MHz (frequency band Europe)
    SpeedLight  = 3e8;         % Speed of light
    eta = 2;                   % expoente de atenuação do caminho
    wavelength = SpeedLight/Freq_Band;

    % definimos as portadoras a serem utilizadas target_pattern
    pattern = zeros(length(TimeStamp),fragment_length);
    pr_h_g_D = zeros(length(TimeStamp),fragment_length);
    pack_tx_segments=zeros(length(TimeStamp),fragment_length);

    for pack=1:1:length(TimeStamp)
        %% cálculo distância entre os dispositivos e satélite
        rho = sqrt(rand(1,1)*(max(Ground_distance)^2));
        Theta = rand(1,1)*2*pi;                      
        Coordinates(:,1) = cos(Theta).*rho;
        Coordinates(:,2) = sin(Theta).*rho;
        distance = sqrt(Coordinates(:,1)^2 + Coordinates(:,2)^2)';
        kD = k(randi(length(k))); % fator Riciano
 
        for frag = 1:1:fragment_length 
            if frag == 1
                pack_tx_segments(pack,frag) = TimeStamp (pack);
            elseif frag > 1 && frag <=(header_size+1)
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + Header_duration;
               
            elseif frag == (header_size+2)
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + Transceiver_wait;

            else
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + fragment_duration;
            end

            % definição das portadoras em cada tx
            if frag == 1
                pattern(pack,frag) = randi(OBW_channels,1,1);
            elseif frag == (header_size+1)
                pattern(pack,frag) = 0; % não exist oprtadora twait
            else
                dif_track=0;
                while dif_track < 8  % 8 x 488 = 3.9 kHz spacing
                    pattern(pack,frag) = randi(OBW_channels,1,1);
                    dif_track=abs(pattern(pack,frag)-pattern(pack,frag-1));
                end
            end

            % potência recebida para cada segmento
            muD = sqrt(kD/(2*(kD+1)));          % Mean 
            sigmaD = sqrt(1/(2*(kD+1)));        % Variance                          
            hrD=sigmaD*randn(1,1)+muD;          % parte real coeficiente de desvanecimento
            hiD=1j*(sigmaD*randn(1,1)+muD);     % aprte imaginária coeficiente de desvanecimento
            h1D=(abs(hrD+hiD)).^2;    % módulo do coeficiente de desvanecimento Riciano
            pr_h_g_D(pack,frag) = Pt*h1D*Gr*Gt*((wavelength/(4*pi*distance))^eta);
        end
    end
end
