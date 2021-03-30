% plot damper
function PlotDamper(damp_x,damp_y,damp_equ)
    % ------------------------------
    % damp_equ: when SMD system is on equilibrium point, the rightest x coordinate of damper
    % damp_origin: origin x coordinate of damper
    % x_range: width of damper
    % y_range: height of damper
    % dlen: the length of right/left horizontal part
    % x_left: x coordinate of damper(left)
    % x_right: x coorindate of damper(right)
    % ver_range: length of short vertical line of damper
    % -------------------------------
    damp_origin = 0;
    x_range = 2.5; y_range = 1.5;
    dlen = 0.5;

    % ----- plot damper -----
    x_left = damp_x - dlen - x_range; 
    x_right = damp_x - dlen;
    y_up = damp_y + y_range/2; y_down = damp_y - y_range/2;
    up = [x_left, x_right, y_up, y_up];
    right = [x_right,x_right,y_up,y_down];
    down = [x_right,x_left,y_down,y_down];
    range = [up;right;down];
    plot([range(:,1),range(:,2)],[range(:,3),range(:,4)],'k','LineWidth',2); hold on;
    % -----------------------
    
    % ----- plot vertical line ------
    ver_range = 1; ver_center = damp_equ - dlen - x_range/2;
    ver_up = damp_y + ver_range/2; ver_down = damp_y - ver_range/2;
    plot([ver_center,ver_center],[ver_down,ver_up],'k','LineWidth',2); hold on;
    % -------------------------------
    
    % ----- plot horizontal line -----
    hor_left = damp_x - dlen;
    plot([damp_origin,ver_center],[damp_y,damp_y],'k','LineWidth',2); hold on;
    plot([hor_left,damp_x],[damp_y,damp_y],'k','LineWidth',2); hold off;
    % --------------------------------
end