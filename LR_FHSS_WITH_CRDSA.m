

function [pattern] = LR_FHSS_WITH_CRDSA(pack_tx_segments, OBW_channels, Header_N)

    %% gerando parâmetros utilizados na simulação
    % definindo as portadoras

    
    %% Frequency-time scheduling of target transmission: which fragment is using the specific channel for a specific time?
    % definimos as portadoras a serem utilizadas target_pattern
    pattern=zeros(1,size(pack_tx_segments,2));

    for row=1:1:size(pack_tx_segments,1)
        for col=1:1:size(pack_tx_segments,2)
            if col == 1
                pattern(row, col) = randi(OBW_channels,1,1);
            else
                dif_track=0;
                while dif_track < 8  % 8 x 488 = 3.9 kHz spacing
                    pattern(row, col) = randi(OBW_channels,1,1);
                    dif_track=abs(pattern(pack,frag)-pattern(pack,frag-1));
                end
            end
        end
    end

end
