% Função para gerar as portadoras, distâncias, potência recebida e o
% coeficiente de desvanecimento para cada transmissão de fragmento

function [x] = Analysys(pattern, pr, pack_tx_segments, Simulation_T, Header_size, Last_fragment_duration)
    t_window = 2*(pack_tx_segments(1, size(pack_tx_segments,2)));
    delta_window = t_window/2;
    Header_duration = 0.233; %233 ms long headers
    F_duration = 0.05;       %50 ms payload data fragments

    % rever este valor, funciona somente para DR8, limiar para decodificar
    % pacote
    Threshold = 4;

    t_window_start = 0;

    % removendo twait (desnecessário para análise)
    pack_tx_segments(:,Header_size+1) = [];
    pattern(:,Header_size+1) = [];
    pr(:,Header_size+1) = [];

    %% separando a janela de tempo detrerminada
    while t_window_start < Simulation_T

        %% selecionando todos os nós transmitidos dentro da janela
        [index_window_row, index_window_col] = find(pack_tx_segments>=t_window_start & pack_tx_segments<=(t_window_start+t_window));
        index_window = [index_window_row, index_window_col];
        keys = unique(index_window_row);
        
        %% capturando todos os fragmentos, portadoras e portências dentro da janela
        pack_tx_segments_window = pack_tx_segments(keys,:);
        pack_tx_segments_window(pack_tx_segments_window <= t_window_start) = 0;
        pack_tx_segments_window(pack_tx_segments_window >= (t_window_start+t_window)) = 0;

        pattern_window = (pack_tx_segments_window ~= 0) .* pattern(keys,:);
        pr_window = (pack_tx_segments_window ~= 0) .* pr(keys,:);

        %% identificando headers na trasmissão
        [hr, hc] = find(pack_tx_segments_window(:, 1:3) ~= 0);
        header_tx_window = unique(hr);

        rf = (pack_tx_segments_window ~= 0);

        %% inicio da análise por pacotes
        for i=1 : length(header_tx_window)
            
            node_tx = pack_tx_segments_window(header_tx_window(i),:);
            % primeira verificação, checa se chegaram headers o suficiente
            % e se existe fragmentos suficiente para decodificar o pacote
            if (sum(node_tx(1:3)~=0)>1 && sum(node_tx(1,4:end)~=0)>Threshold)

                % criando cópia das tx dentro da janela para eliminar o
                % dispositivo transmitido e iniciar a localização das
                % colisões
                pack_tx = pack_tx_segments_window;
                pack_tx(header_tx_window(i),:) = -10;
                pack_tx(pack_tx==0) = -10;
                transmission_sim = [];

                % início de análise para colisões
                for seg=1:1:length(node_tx)

                    % colisões durante transmissão header
                    if(seg ~= length(node_tx) && seg <= Header_size)
                        [header_sim_row, header_sim_col] = find(pack_tx(:,1:Header_size) >= (node_tx(seg)-Header_duration) ...
                            & pack_tx(:,1:Header_size) <= (node_tx(seg)+Header_duration));

                        [frag_sim_row, frag_sim_col] = find(pack_tx(:,(Header_size+1):end-1) >= (node_tx(seg)-F_duration) ...
                            & pack_tx(:,(Header_size+1):end-1) <= (node_tx(seg)+Header_duration));

                        [last_sim_row, last_sim_col] = find(pack_tx(:,end) >= (node_tx(seg)-Last_fragment_duration) & ...
                            pack_tx(:,end) <= (node_tx(seg)+Header_duration));
                    
                    % colisões durante transmissão fragmento
                    elseif(seg~=length(node_tx) && seg > Header_size)
                        [header_sim_row, header_sim_col] = find(pack_tx(:,1:Header_size) >= (node_tx(seg)-Header_duration) ...
                            & pack_tx(:,1:Header_size) <= (node_tx(seg)+F_duration));
                        
                        [frag_sim_row, frag_sim_col] = find(pack_tx(:,(Header_size+1):end-1) >= (node_tx(seg)-F_duration) ...
                            & pack_tx(:,(Header_size+1):end-1) <= (node_tx(seg)+F_duration));

                        [last_sim_row, last_sim_col] = find(pack_tx(:,end) >= (node_tx(seg)-Last_fragment_duration) ...
                            & pack_tx(:,end) <= (node_tx(seg)+F_duration));

                    % colisões durante o último fragmento
                    else
                        [header_sim_row, header_sim_col] = find(pack_tx(:,1:Header_size) >= (node_tx(seg)-Header_duration) ...
                            & pack_tx(:,1:Header_size) <= (node_tx(seg)+Last_fragment_duration));
                        
                        [frag_sim_row, frag_sim_col] = find(pack_tx(:,(Header_size+1):end-1)>=(node_tx(seg)-F_duration) ...
                            & pack_tx(:,(Header_size+1):end-1)<=(node_tx(seg)+Last_fragment_duration));
                        
                        [last_sim_row, last_sim_col] = find(pack_tx(:,end) >= node_tx(seg) ...
                            & pack_tx(:,end)<=(node_tx(seg)+Last_fragment_duration));
                
                    end

                    segmento_header = [];
                    segmento_frag = [];
                    segmento_last = [];

                    segmento_header(1:length(header_sim_row),1) = seg;
                    segmento_frag(1:length(frag_sim_row),1) = seg;
                    segmento_last(1:length(last_sim_row),1) = seg;

                    transmission_sim = [transmission_sim; header_sim_row, header_sim_col, segmento_header; 
                        frag_sim_row, frag_sim_col, segmento_frag; 
                        last_sim_row, last_sim_col, segmento_last];
                
                end

                transmission_sim = unique(transmission_sim,'rows')


                
            end
            
        end


        % incrementando delta para formação de nova janela
        t_window_start = t_window_start + delta_window;
    end
    
    x = 0;
end





























