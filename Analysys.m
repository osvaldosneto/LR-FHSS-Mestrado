% Função para gerar as portadoras, distâncias, potência recebida e o
% coeficiente de desvanecimento para cada transmissão de fragmento

function [x] = Analysys(pattern, distance, pr_h_g_D, h1D, pack_tx_segments, Simulation_T, Header_size)
    t_window = 2*(pack_tx_segments(1, size(pack_tx_segments,2)));
    delta_window = t_window/2;

    t_window_start = 0;

    % removendo twait (desnecessário para análise)
    pack_tx_segments(:,Header_size+1) = [];
    pattern(:,Header_size+1) = [];
    pr_h_g_D(:,Header_size+1) = [];
    h1D(:,Header_size+1) = [];

    while t_window_start < Simulation_T
        
        
        % selecionando todos os nós transmitidos dentro da janela
        [index_window_row, index_window_col] = find(pack_tx_segments>=t_window_start & pack_tx_segments<=(t_window_start+t_window));
        index_window = [index_window_row, index_window_col];
        
        % capturando todos os fragmentos dentro da janela
        keys = unique(index_window_row);
        pack_tx_segments_window = zeros(size(keys,1), size(pack_tx_segments,2));
        for i = 1:numel(keys)
            group_index = (index_window_row == keys(i));
            pack = pack_tx_segments(keys(i),index_window_col(group_index));
            
            % ajustando size do vetor para adicionar na janela
            pack = padarray(pack, [0, size(pack_tx_segments_window, 2) - size(pack, 2)], 'post');
            
            pack_tx_segments_window(i, :) = pack;
        end


        % incrementando delta para formação de nova janela
        t_window_start = t_window_start + delta_window;
    end
    
    x = 0;
end