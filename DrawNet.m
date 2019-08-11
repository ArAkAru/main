function net = DrawNet(net,N,CH,D,SX,SY)

% normal nodes
%figure();
%hold on;
plot(net(2,~CH&~D),net(3,~CH&~D),'ko','MarkerSize',...
    5,'MarkerFaceColor','b'); hold on;
% dead nodes
plot(net(2,D),net(3,D),'ko','MarkerSize',...
    5,'MarkerFaceColor','k'); hold on;
% Cluster heads
plot(net(2,CH),net(3,CH),'ko','MarkerSize',...
    5,'MarkerFaceColor','r'); hold on;
% The sink
scatter(SX,SY,100,'mo','filled')
text(SX+2,SY+2,'Base','FontSize',10,'VerticalAlignment','Baseline');
% titles
s = int2str((1:N)');
text(net(2,:)+1,net(3,:)+1,s,'FontSize',8,'VerticalAlignment','Baseline');
xlabel('\it x \rm [m] \rightarrow');
ylabel('\it y \rm [m] \rightarrow');

hold off