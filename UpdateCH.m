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
%���������� ��� �� ������ �������� 1 ��� -� ����� �� ����������� � ��� ����
%������ ������� ����� ��������� ���,
T(idx) = p / (1-p * mod(r, round(1 / p)));%��� ������ ��� ����� ������
T(~idx) = 0;%���� ��� ��� �������,������ �� ����� 
T(D) = 0;%���� ���� �� ����� �������
CH = tmp<T; %�� ���� ��������� ����� ������� ������� �����

G(CH) = 0; 
