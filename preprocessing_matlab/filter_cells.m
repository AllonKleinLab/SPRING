function [Efiltered,cell_filter] = filter_cells(E,min_reads)
% Filter out cells < min_reads UMIs
%
% Inputs
%   E           = Expression matrix. Rows are cells and columns are genes.
%   min_reads   = Minimum number of UMIs to keep a cell
% 
% Outputs
%   Efiltered   = Cell-filtered expression matrix
%   cell_filter = Boolean array for filtering cells
%
%%
    total_counts = sum(E,2);
	cell_filter = total_counts >= min_reads;
	if sum(cell_filter) == 0
        Efiltered = [];
    else
        Efiltered = E(cell_filter,:);
    end
end