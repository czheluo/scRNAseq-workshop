colo(1,1)='b';colo(1,2)='b';
colo(2,1)='g';colo(2,2)='g';
colo(3,1)='c';colo(3,2)='c';
colo(4,1)='r';colo(4,2)='r';
colo(5,1)='m';colo(5,2)='m';
colo(6,1)='y';colo(6,2)='y';
colo(7,1)='k';colo(7,2)='k';
colo(8,1)='w';colo(8,2)='b';
colo(9,1)='w';colo(9,2)='g';
colo(10,1)='w';colo(10,2)='c';
colo(11,1)='w';colo(11,2)='r';
colo(12,1)='w';colo(12,2)='m';
colo(13,1)='w';colo(13,2)='k';
colo(14,1)='k';colo(14,2)='g';
colo(15,1)='k';colo(15,2)='c';
colo(16,1)='k';colo(16,2)='r';
colo(17,1)='k';colo(17,2)='m';
colo(18,1)='r';colo(18,2)='b';
colo(19,1)='r';colo(19,2)='g';
colo(20,1)='r';colo(20,2)='c';
colo(21,1)='r';colo(21,2)='k';
colo(22,1)='k';colo(22,2)='c';
colo(23,1)='k';colo(23,2)='m';
colo(24,1)='m';colo(24,2)='g';
colo(25,1)='c';colo(25,2)='b';
mk(1:7)='o';mk(8:14)='s';mk(15:21)='p';mk(22:25)='^';
if k==1
    hold off
    plot(x1(:,1),.95,'o')
    hold on
    for i=1:g
        for l=1:n
            if sum(abs(ym(l,:,i)))>0
                plot(x1(l,1),i*.8/g,mk(i),'markerfacecolor',colo(i,1),'markeredgecolor',colo(i,2),'linewidth',2.3)
                axis([x1a,x1b,0,1])
                %text(x1(l,1),i*.8/g+(rand(1,1)-.45)*.1,num2str(l),'fontsize',10)
            end
        end
    end
elseif k==2
    hold off
    plot(x1(:,1),x1(:,2),'o')
    hold on
    for i=1:g
        for l=1:n
            if sum(abs(ym(l,:,i)))>0
               plot(x1(l,1),x1(l,2),mk(i),'markerfacecolor',colo(i,1),'markeredgecolor',colo(i,2),'linewidth',2.3)
               %text(x1(l,1),x1(l,2)+.015*r3,num2str(l),'fontsize',11,'fontweight','b')
            end
        end
    end
elseif k==3
    hold off
    plot3(x1(:,1),x1(:,2),x1(:,3),'.')
    hold on
    for i=1:g
        for l=1:n
            if sum(abs(ym(l,:,i)))>0
                plot3(x1(l,1),x1(l,2),x1(l,3),mk(i),'markerfacecolor',colo(i,1),'markeredgecolor',colo(i,2),'linewidth',2.5)
                %text(x1(l,1),x1(l,2),x1(l,3)+.015*r3,num2str(l),'fontsize',10)
            end
        end
        %stem3(x1(:,1),x1(:,2),x1(:,3))
        scatter3(x1(:,1),x1(:,2),x1(:,3))
    end
    grid on
    xlabel('PC1'),ylabel('PC2'),zlabel('PC3')
    axis tight
end
