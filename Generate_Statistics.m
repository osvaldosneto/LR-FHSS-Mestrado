function [d_plot, prob] = Generate_Statistics(sim, distance)

    %% extraindo estatísticas da transmissão
    d_slice = [1.1e6 1.3e6 1.5e6 1.7e6 1.9e6 2.1e6 2.3e6];
    dif1 = 0; devices1 = 0;
    dif2 = 0; devices2 = 0;
    dif3 = 0; devices3 = 0;
    dif4 = 0; devices4 = 0;
    dif5 = 0; devices5 = 0;
    dif6 = 0; devices6 = 0;
    dif7 = 0; devices7 = 0;
    
    for n=1:1:size(sim,1)
        for d=1:1:length(d_slice)
    
            switch d
                case 1
                    devices = find(distance < d_slice(d));
                    dif1 = dif1 + length(devices)-length(setdiff(devices,sim(n,:)));
                    devices1 = devices1 + length(devices);
    
                case 2
                    devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                    dif2 = dif2 + length(devices)-length(setdiff(devices,sim(n,:)));
                    devices2 = devices2 + length(devices);
    
                case 3
                    devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                    dif3 = dif3 + length(devices)-length(setdiff(devices,sim(n,:)));
                    devices3 = devices3 + length(devices);
    
                case 4
                    devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                    dif4 = dif4 + length(devices)-length(setdiff(devices,sim(n,:)));
                    devices4 = devices4 + length(devices);
    
                case 5
                    devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                    dif5 = dif5 + length(devices)-length(setdiff(devices,sim(n,:)));
                    devices5 = devices5 + length(devices);
    
                case 6
                    devices = find(distance >= d_slice(d-1) & distance < d_slice(d));
                    dif6 = dif6 + length(devices)-length(setdiff(devices,sim(n,:)));
                    devices6 = devices6 + length(devices);
    
                case 7
                    devices = find(distance >= d_slice(d));
                    dif7 = dif7 + length(devices)-length(setdiff(devices,sim(n,:)));
                    devices7 = devices7 + length(devices);
            end
    
        end
    
    end
    
    dif1 = dif1/devices1;
    dif2 = dif2/devices2;
    dif3 = dif3/devices3;
    dif4 = dif4/devices4;
    dif5 = dif5/devices5;
    dif6 = dif6/devices6;
    dif7 = dif7/devices7;
    
    d_plot = [1e6 1.2e6 1.4e6 1.6e6 1.8e6 2.0e6 2.2e6];
    prob = [dif1 dif2 dif3 dif4 dif5 dif6 dif7];

end