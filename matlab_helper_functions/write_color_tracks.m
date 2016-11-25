function write_color_tracks(ctracks, fname)
% Write color data to csv file
%
% ctracks  = Cell array. Each row is a color track. The first entry in the
%            row is the name of the track and subsequent entries give 
%            specific color values. If there are N cells, then this cell
%            array should have N+1 columns. 
% fname    = Path to save the color tracks file. 
%
%%
    for i=1:size(ctracks,1)
        for j=2:size(ctracks,2)
            val = ctracks{i,j};
            val = round(val,1);
            val = num2str(val);
            val = strrep(val,'0.','.');
            ctracks{i,j} = val;
        end
    end
	cell2csv(fname,ctracks);
end