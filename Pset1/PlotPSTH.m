% 1.2

function PlotPSTH(cell_array, bin_num)

    edges = linspace(0, .6, bin_num); 
    psth = zeros(1,length(edges));

    for j = 1:length(cell_array)
        psth = psth + histc(cell_array{j}, edges);
    end

    avg = psth / length(cell_array);
    plot(edges, avg);
    xlim([0 .6]);
    xlabel('Time (sec)');
    ylabel('# of spikes');
end
