
function [Efiltered,cell_filter] = filter_cells(E,min_reads)
    total_counts = sum(E,2);
	cell_filter = total_counts >= min_reads;
	if sum(cell_filter) == 0
        Efiltered = [];
    else
        Efiltered = E(cell_filter,:);
    end
end