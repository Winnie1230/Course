function PlotMap(map_pos,map_len)
    for i = 1:length(map_pos)
        rectangle('Position',[map_pos(i,:) map_len(i,:)])
    end
end