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
        
        index_window = sortrows(index_window, 1); % index de todos os nós transmitindo na mesma janela de tempo
        

        t_window_start = t_window_start + delta_window;
    end
    
    x = 0;
end