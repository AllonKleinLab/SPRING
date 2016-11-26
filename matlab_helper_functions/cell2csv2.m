function cell2csv2(file,cellArray)
%file : string with the name of the file to create
%cellArray : cell array of numbers and strings
delimiter = ',';
types = cellfun(@class,cellArray,'UniformOutput',false);
inds = strcmp(types,'char');
types(inds)={'%s'};
types(~inds)={'%f'};
types = strcat(types,delimiter);
%Write line by line
fid = fopen(file,'w');
for indLine = 1:size(cellArray,1)
formatLine = strcat(types{indLine,:});
formatLine(end) = [];%delete the trailing ;
fprintf(fid,[formatLine,'\n'],cellArray{indLine,:});
end
fclose(fid);