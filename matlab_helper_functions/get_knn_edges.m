function edges = get_knn_edges(dmat, k)
	edges = [];
    for i = 1:size(dmat,1)
    	[dSorted dIdx] = sort(dmat(i,:));
        for j = dIdx(1:k)
            edge = sort([i,j]);
            if size(edges,1) == 0 || sum(ismember(edges, edge, 'rows'))==0
                edges = [edges; edge];
            end
        end
    end
end