function cellarray = struct_field_qualified(cellarray)
    for i = 1:size(cellarray,1)
        for j = 1:size(cellarray,2)
            if ischar(cellarray{i,j})
                if ~isletter(cellarray{i,j}(1))
                    cellarray{i,j} = strcat('A',cellarray{i});
                end
                cellarray{i,j} = strrep(cellarray{i,j},'-','_');
                cellarray{i,j} = strrep(cellarray{i,j},'.','_');
                cellarray{i,j} = strrep(cellarray{i,j},'/','_');
                cellarray{i,j} = strrep(cellarray{i,j},' ','_');
            end
        end
    end
end