% Função para gerar as portadoras, distâncias, potência recebida e o
% coeficiente de desvanecimento para cada transmissão de fragmento

function [x] = Analysys(pattern, pr, pack_tx_segments, Simulation_T, Header_size)
    t_window = 2*(pack_tx_segments(1, size(pack_tx_segments,2)));
    delta_window = t_window/2;
    % rever este valor, funciona somente para DR8
    Threshold = 4;

    t_window_start = 0;

    % removendo twait (desnecessário para análise)
    pack_tx_segments(:,Header_size+1) = [];
    pattern(:,Header_size+1) = [];
    pr(:,Header_size+1) = [];

    %% separando a janela de tempo detrerminada
    while t_window_start < Simulation_T

        % selecionando todos os nós transmitidos dentro da janela
        [index_window_row, index_window_col] = find(pack_tx_segments>=t_window_start & pack_tx_segments<=(t_window_start+t_window));
        index_window = [index_window_row, index_window_col];
        keys = unique(index_window_row);
        
        % capturando todos os fragmentos dentro da janela
        pack_tx_segments_window = pack_tx_segments(keys,:);
        pack_tx_segments_window(pack_tx_segments_window <= t_window_start) = 0;
        pack_tx_segments_window(pack_tx_segments_window >= (t_window_start+t_window)) = 0;

        
        % identificando headers na trasmissão
        headers_frag = index_window(index_window(:,2) <= 3);
        

        % incrementando delta para formação de nova janela
        t_window_start = t_window_start + delta_window;
    end
    
    x = 0;
end