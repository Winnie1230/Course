function PlotCar(center, dir)
    car_w = 20; % cm
    car_l = 15; % cm
    car_dir = [0;1];

    corner = [center(1)-car_w/2, center(2)-car_l/2; ...
              center(1)+car_w/2, center(2)-car_l/2; ...
              center(1)+car_w/2, center(2)+car_l/2; ...
              center(1)-car_w/2, center(2)+car_l/2];

    R = [cosd(dir) -sind(dir); sind(dir) cosd(dir)]; % rotation matrix
    center = R*center;
    corner = (R*corner')';
    car_dir = R*car_dir;
    factor = 4;
    plot(center(1),center(2),'r+'); hold on
%     plot(sensor(:,1),sensor(:,2),'bo','MarkerSize',15,'MarkerFace','b');
%     plot(corner(:,1), corner(:,2),'bo','MarkerSize',10); hold on
    plot([corner(1,1) corner(2,1)],[corner(1,2) corner(2,2)],'LineWidth',1.5,'Color','k'); hold on
    plot([corner(2,1) corner(3,1)],[corner(2,2) corner(3,2)],'LineWidth',1.5,'Color','k'); hold on
    plot([corner(3,1) corner(4,1)],[corner(3,2) corner(4,2)],'LineWidth',1.5,'Color','k'); hold on
    plot([corner(1,1) corner(4,1)],[corner(1,2) corner(4,2)],'LineWidth',1.5,'Color','k'); hold on

    quiver(center(1),center(2),car_dir(1)*factor,car_dir(2)*factor,0,'LineWidth',2,'Color','g','MaxHeadSize',3);
    axis equal
end