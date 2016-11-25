function edges = get_knn_edges(dmat, k)
% Get k-nearest-neighbor graph edges from a distance matrix
%
% dmat     = Distance matrix. Entry i,j is the distance between nodes i and j
% k        = Number of edges assigned to each node
%
%%
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