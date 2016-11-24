function [Efiltered, gene_filter] = filter_genes(E, Ecutoff, Vcutoff)
	mean_filter = mean(E,1) >= Ecutoff;
	var_filter = var(E,1) ./ (mean(E,1)+.0001) >= Vcutoff;
	gene_filter = all([mean_filter; var_filter],1);
	Efiltered = E(:,gene_filter);
end