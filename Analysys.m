% Função para gerar as portadoras, distâncias, potência recebida e o
% coeficiente de desvanecimento para cada transmissão de fragmento

function [receive_success, colision] = Analysys(pattern, pr, pack_tx_segments, Header_size, Last_fragment_duration, SIR_op)
    t_window = 2*(pack_tx_segments(1,end)-pack_tx_segments(1,1));
    delta_window = t_window/2;
    Header_duration = 0.233; %233 ms long headers
    F_duration = 0.05;       %50 ms payload data fragments

    % rever este valor, funciona somente para DR8, limiar para decodificar
    % pacote
    Threshold = 4;
    t_window_start = 0;
    t_packtx_segments_end = sort(pack_tx_segments);
    t_end = t_packtx_segments_end(end);

    % armazenamento dos pacotes recebidos por completo
    receive_success = [];
    colision = 0;

    % removendo twait (desnecessário para análise)
    pack_tx_segments(:,Header_size+1) = [];
    pattern(:,Header_size+1) = [];
    pr(:,Header_size+1) = [];

    %% janela de tempo determinada
    while t_window_start < t_end

        % capturando todos os fragmentos, portadoras e portências dentro da janela
        pack_tx_segments_window = pack_tx_segments;
        pack_tx_segments_window(pack_tx_segments_window <= t_window_start) = 0;
        pack_tx_segments_window(pack_tx_segments_window >= (t_window_start+t_window)) = 0;
    
        pattern_window = (pack_tx_segments_window ~= 0) .* pattern;
        pr_devices = (pack_tx_segments_window ~= 0) .* pr;
        
        
        % validando fim de loop
        loop = 1;
        receive_success_history = [];
        
        while loop == 1

            pr_I = zeros(size(pr_devices));
        
            %% identificando headers na trasmissão
            [hr, hc] = find(pack_tx_segments_window(:, 1:3) ~= 0);
            header_tx_window = unique(hr);
            header_tx_window = setdiff(header_tx_window,receive_success);
        
            %% inicio da análise por pacotes tendo headers como fonte de pesquisa
            for i=1 : length(header_tx_window)
                    
                node_tx = pack_tx_segments_window(header_tx_window(i),:);
        
                % criando cópia das tx dentro da janela para eliminar o
                % dispositivo transmitido e iniciar a localização das
                % transmissões simultâneas
                pack_tx = pack_tx_segments_window;
                pack_tx(header_tx_window(i),:) = -10;
                pack_tx(pack_tx==0) = -10;
        
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
        
                    transmission_sim = [header_sim_row, header_sim_col; 
                             frag_sim_row, frag_sim_col + Header_size; %soma para ajuste do index devido a função find eliminar header
                             last_sim_row, last_sim_col + size(pack_tx_segments_window,2)-1]; %soma para ajuste do index
        
                    % verificando colisão entre portadoras
                    for c=1:1:size(transmission_sim(:,1))
        
                        % verifica se houve colisão
                        if pattern_window(transmission_sim(c,1), transmission_sim(c,2)) == pattern_window(header_tx_window(i), seg)
                            pr_I(header_tx_window(i), seg) = pr_I(header_tx_window(i), seg) + pr_devices(transmission_sim(c,1), transmission_sim(c,2));
                            colision = colision+1;
                        end
        
                    end
        
                end
              
            end
        
            %% análise de pacotes recebidos
            pack_collided_decoded = (pr_I ~= 0);    % armazenamento da posição onde houve a colisão
            pr_devices_collided = pr_devices .* pack_collided_decoded;
        
            % relação Sinal interferência
            sir_pack_received = pr_devices_collided./abs(pr_I-pr_devices_collided);
            sir_pack_received = fillmissing(sir_pack_received, 'constant', 0);
        
            % identificando pacotes recebidos 
            % 1 recebidos com sucesso 
            % 0 caso contrário
            received_pack = (sir_pack_received>4 | sir_pack_received==0) .* (pack_tx_segments_window ~= 0); % multiplicando pela janela
                
            % selecionando pacotes recebidos com sucesso, verifica se recebeu
            % ao menos um header e ao menos 4 fragmentos
            pack_success = find(sum(received_pack(:,1:3),2) > 0 & sum(received_pack(:,4:end),2) >= Threshold);
        
            receive_success = [receive_success; pack_success];
            receive_success = unique(receive_success);
               
            if SIR_op
                if isequal(receive_success, receive_success_history)
                    loop = 0;
                end
                % removendo recebidos da análise
                pr(receive_success,:) = 0;
                receive_success_history = receive_success;
            else
                loop = 0;
            end

        end

        % incrementando delta para formação de nova janela
        t_window_start = t_window_start + delta_window;

    end

end
