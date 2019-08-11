function d = fun(x,y,net,CH)
%disp(net(2,CH));
if isempty(net(2,CH))%если головных нету.дистанция 0
    %disp('s');
    d=0;
else
    d=min(sqrt((net(2,CH) - x).^2 + (net(3,CH) - y).^2));%минимальное расстояние между головным и обычным узлом
end