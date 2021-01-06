% plot SMD system by inputing the mass center position
function PlotSMD(mass_center,equ_pt,mass_width,mass_height,x_range)
    x_axis = 7/4*(x_range+1) + x_range;
    y_axis = 6;
    % ----- plot Wall -----
    plot([0 x_axis],[0 0],'k','LineWidth',2); hold on;
    plot([0 0],[0 y_axis],'k','LineWidth',2); hold on;

    PlotMass(mass_center,mass_width,mass_height); hold on;
    % ----- plot spring -----
    spr_x = mass_center(1)-mass_width/2; spr_y = 3.5;
    PlotSpring(spr_x,spr_y,x_range/3); hold on;
    
    % ----- plot damper -----
    damp_x = spr_x; damp_y = 1;
    damp_equ = equ_pt(1)-mass_width/2;
    PlotDamper(damp_x,damp_y,x_range,damp_equ); hold off;
    axis([0 x_axis 0 y_axis]);
    axis equal tight;
end