function PlotQuiver(p1,dp,c)
%     dp = p2 - p1
    quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'Color',c,'MaxHeadSize',3);
end