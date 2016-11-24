function write_color_tracks(ctracks, fname)
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