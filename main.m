% Clean memory and command window
clear,clc,close all

%% Parameters
N = 50;          %Number of nodes
W = 500;          %length of the network
L = 500;          %width of the network
p = 0.1;          %desired percentage of cluster heads
R = 10;           %radius of radio range for a node
num_rounds = 200; %Max Number of simulated rounds

Ei = 0.1;              % Initial energy of each node
Etrans = 70.0000e-09; % Energy for transmitting one bit 
Erec = 70.0000e-09;   % Energy for receiving one bit 
Eagg = 5.0000e-09;   % Data aggregation energy
Efs = 10.0000e-12;
Emp=0.0013e-12;
d0=sqrt(Efs/Emp);
CHpl = 4096;         % Packet size for cluster head per round (bits)
NonCHpl = 1500;       % Packet size for normal node per round (bits)

Tsetup = 3;          % average Time in seconds taken in setup phase
Tss = 10;            % average Time in seconds taken in steady state phase

% Position of sink
SX = 250; SY = 250;
% Preallocation for energy calculations
E = Ei*ones(1,N);          % Energy left in each node
EH = zeros(1,num_rounds); 
% Preallocation for dead nodes calculations
Ecrit = 0;                 % Critical energy left in node to call it alive
Ncrit = fix((95/100)*N);   % Critical number for dead nodes to stop simulation
Dead = false(1,N);         % Status of all nodes 0:Alive 1:Dead
DeadH = zeros(1,num_rounds);
% Preallocation for Bits sent calculations
BitsH = zeros(1,num_rounds);


%% Building the WSN
% 1st row: states of being a CH, 1:never been CH, 0:has been CH 
% 2nd: x-position, 3rd: y-position
net = [ones(1,N);rand([1,N])*W;rand([1,N])*L];
figure('Position',[34 30 792 613]);

%% Simulating LEACH for each round
for r=1:num_rounds % iterating on each round
    CH=1<0;
    %%%% Choosing Clusters heads %%%%
    %[net(1,:),CH] = UpdateCH(net(1,:),Dead,p,r);
   % disp(CH);
    %%%% Energy calculations %%%%
    EH(r) = sum(E); %get all energy left in all nodes
    % first CH
    %numClust = length(find(CH));%кол-во головных узлов
    %D = sqrt((net(2,CH) - SX).^2 + (net(3,CH) - SY).^2);%расстояние меежду головным и бс
    %E(CH) = E(CH) - (((Etrans+Eagg)*CHpl)+(Efs*CHpl*(D.^ 2))+(NonCHpl*Erec*round(N/numClust)));
    % second rest of nodes
    rest = N-sum(double(Dead));
    %mD = zeros(1,rest); 
    tmp = net(2:3,~Dead);%берем координаты узлов которые не были головными и не умерли
    m=1;
    MasIDX=zeros(1,N);
    
    D =  sqrt((tmp(1,:) - SX).^2 + (tmp(2,:) - SY).^2); %высчитали минимальное расстояние между головным и обычнчм узлом
 
    
    for j=1:length(D)
    if(D(j)<d0)
    E(j)=E(j)-((Etrans*NonCHpl)+(Efs*NonCHpl*(D(j)^ 2)));
    end
    if(D(j)>=d0)
        E(j)=E(j)-((Etrans*NonCHpl)+(Emp*NonCHpl*(D(j)^ 4)));
    end
    end
    
    
     %E(~Dead&dist) = E(~Dead&dist)-((Etrans*NonCHpl)+(Efs*NonCHpl*(Def.^ 2)));
    
        
    
   % E(~Dead&~dist) = E(~Dead&~dist)-((Etrans*NonCHpl)+(Emp*NonCHpl*(Dmp.^ 4)));
    %k=(Efs*NonCHpl*u);
    
    %E(~Dead) = E(~Dead) - ((NonCHpl*Etrans)*mD + ((Erec+Eagg)*CHpl));
    %finally updating alive status to all nodes
    E(Dead) = 0;
    Dead(E<=Ecrit) = true ; DeadH(r)=sum(double(Dead));
    %%%% sent bits %%%%
    BitsH(r+1) =  BitsH(r)+rest*NonCHpl;
    
    %%%% Showing updated net %%%%
    net = DrawNet(net,N,CH,Dead,SX,SY);
    title(['Normal nodes:Black ---- CH:Red ---- Dead:Empty circle --- round (',num2str(r),')']);
    drawnow
    
    if DeadH(r)>=Ncrit,break;end % Stop simulation when 5% or less is alive
end

%% Plotting analysis of network preformance

%% Plotting analysis of network preformance

T = (Tsetup+Tss)*(0:r-1);
EH = EH(1:r); EHdis = EH;
DeadH = DeadH(1:r); AliveH = N-DeadH;
BitsH = BitsH(2:r+1);

figure('Position',[131 59 792 613]);
plot(T,EHdis,'-x'); xlabel('Time (s)'); ylabel('Energy (j)')
title('Total energy dissipated')
legend('LEACH')

figure('Position',[298 66 792 613]);
plot(T,AliveH,'-x'); xlabel('Time (s)'); ylabel('No of live Sensors nodes')
title('Life time of sensor nodes')
legend('LEACH')

figure('Position',[511 61 792 613]);
plot(T,BitsH,'-x'); xlabel('Time (s)'); ylabel('Throughput (no of packets)')
title(['Throughput (' num2str(N) ' Sensor nodes)'])
legend('LEACH')



