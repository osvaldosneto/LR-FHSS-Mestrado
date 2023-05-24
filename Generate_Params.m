% Função para gerar as portadoras, distâncias, potência recebida e o
% coeficiente de desvanecimento para cada transmissão de fragmento

function [pattern, distance, pr_h_g_D, h1D] = Generate_Params(pack_tx_segments, ...
    OBW_channels, Ground_distance, k)

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
    pattern = zeros(size(pack_tx_segments));
    distance = zeros(size(pack_tx_segments,1),1)';
    pr_h_g_D = zeros(size(pack_tx_segments));
    h1D = zeros(size(pack_tx_segments));

    for row=1:1:size(pack_tx_segments,1)
        %% cálculo distância entre os dispositivos e satélite
        rho = sqrt(rand(1,1)*(max(Ground_distance)^2));
        Theta = rand(1,1)*2*pi;                      
        Coordinates(:,1) = cos(Theta).*rho;
        Coordinates(:,2) = sin(Theta).*rho;
        distance(row) = sqrt(Coordinates(:,1)^2 + Coordinates(:,2)^2)';

        %% Cálculo do sinal recebido (valor desejado)
        kD = k(3); % validar se este fator Riciano está correto 
        
        % Rician fading
        for col=1:1:size(pack_tx_segments,2)  
            muD = sqrt(kD/(2*(kD+1)));          % Mean 
            sigmaD = sqrt(1/(2*(kD+1)));        % Variance                          
            hrD=sigmaD*randn(1,1)+muD;          % parte real coeficiente de desvanecimento
            hiD=1j*(sigmaD*randn(1,1)+muD);     % aprte imaginária coeficiente de desvanecimento
            h1D(row, col)=(abs(hrD+hiD)).^2;    % módulo do coeficiente de desvanecimento Riciano
            pr_h_g_D(row, col) = Pt*h1D(row, col)*Gr*Gt*((wavelength/(4*pi*distance(row)))^eta);
        end

        %% Cálculo das portadoras
        for col=1:1:size(pack_tx_segments,2)
            if col == 1
                pattern(row, col) = randi(OBW_channels,1,1);
            else
                dif_track=0;
                while dif_track < 8  % 8 x 488 = 3.9 kHz spacing
                    pattern(row, col) = randi(OBW_channels,1,1);
                    dif_track=abs(pattern(row, col)-pattern(row, col-1));
                end
            end
        end
    end

end
