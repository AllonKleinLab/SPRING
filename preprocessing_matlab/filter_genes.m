function [Efiltered, gene_filter] = filter_genes(E, Ecutoff, Vcutoff)
% Filter genes based on expression level and variability. Only include
% genes with mean > Ecutoff and fano factor > Vcutoff
%
% Inputs
%   E           = Expression matrix. Rows are cells and columns are genes.
%   Ecutoff     = Minimum mean expression
%   Vcutoff     = Minimum fano factor
% 
% Outputs
%   Efiltered   = Gene-filtered expression matrix
%   gene_filter = Boolean array for filtering genes
%
%%
    mean_filter = mean(E,1) >= Ecutoff;
	var_filter = var(E,1) ./ (mean(E,1)+.0001) >= Vcutoff;
	gene_filter = all([mean_filter; var_filter],1);
	Efiltered = E(:,gene_filter);
end