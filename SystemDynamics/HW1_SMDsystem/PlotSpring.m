% plot spring
function PlotSpring(spr_x,spr_y)
    % ----------------------------------------------------
    % spr_x: spring rightest x
    % spr_y: spring y coordinate
    % slen: length of horizontal line on the left and right
    % saw_num: number of saw
    % spr_origin: origin x coordinate of spring
    % x_left: spring left x coordinate(saw part)
    % x_right: spring right x coordinate(saw part)
    % y_range: vertical length of spring
    % ----------------------------------------------------
    spr_origin = 0;
    slen = 0.5;  saw_num = 6;
    x_left = spr_origin + slen; x_right = spr_x - slen;
    y_range = 2;
    y_down = spr_y-y_range/2; y_up = spr_y+y_range/2;
    saw_x = linspace(x_left, x_right, saw_num*4+1);
    yy = [y_down,spr_y,y_up,spr_y];

    for i = 1:size(saw_x,2)-1
        range(i,:) = [saw_x(i), saw_x(i+1),yy(mod(i,4)+1),yy(mod(i+1,4)+1)];
    end
    plot([spr_origin, spr_origin+slen],[spr_y, spr_y],'k','LineWidth',2); hold on;
    plot([range(:,1),range(:,2)],[range(:,3),range(:,4)],'k','LineWidth',2); hold on;
    plot([x_right, x_right+slen],[spr_y, spr_y],'k','LineWidth',2); hold off;
end