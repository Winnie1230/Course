% plot mass
function PlotMass(mass_center,width,height)
    % width: width of mass
    % height: height of mass
%     width = 2; height = 5;
    plot(mass_center(1),mass_center(2),'r+'); hold on;
    rectangle('Position',[mass_center(1)-width/2 mass_center(2)-height/2 width height],'LineWidth',2);
end