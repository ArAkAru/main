function [G,CH] = UpdateCH(G,D,p,r)

tmp = rand(size(G)); %random number for each node

T = G; 
%disp(T);
idx = T==1;
%disp(idx);
if(mod(r, round(1 / p))==0)%if isempty(find(idx, 1))
    %disp('as');
    G=ones(size(G));
    idx = G==1;
end
%переделать так не пойдет остается 1 нод -в итоге не обнулсяется и его шанс
%именно ппройти порог нереально мал,
T(idx) = p / (1-p * mod(r, round(1 / p)));%чем дальше тем порог больше
T(~idx) = 0;%если уже был главным,сейчас не будет 
T(D) = 0;%если умер не будет главным
CH = tmp<T; %из всех возможных стать главным смотрим порог

G(CH) = 0; 
