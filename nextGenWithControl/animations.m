
if rem(k,10) == 0

if animation == 1
    
    

    plot(x/chord,y/chord,'*','LineWidth',1);
    %plot(CPX,CP,'black','LineWidth',1);
    hold on;
    quiver(CPX/chord,CPY/chord,cprx,cpry,1,'black');
    
    hold on;
    xlim([-1 1]);
    hold on
    ylim([-0.5 0.5]);
    hold on
    daspect([1 1 1]);
    
    %hold on;
    %scatter(XM/chord,YM/chord,10,CV,'filled');
    hold on;
    quiver(0,0,CD(k),CL(k),1,'blue','filled','Linewidth',2);
    %hold on;
    %quiver(0,0,0,CL(k),1,'blue','filled');
    hold off;
    
    
    %pbaspect([2 1 1]);
    pause(0.01);
    %scatter(x,y,10);
    %plot(XZ,YZ,'o');
    
end

end