clc;
clear all;
close all;
% Initial state x(0)
X0=[-5;0;(0*pi)/180];
vk=0;
wk=0;
Ts=0.005;
thetak=(0*pi)/180;
D=[-0.01;-0.01;-0.01];
% D=[0;0;0];
N=5;
s=10;
Simlength=(s/Ts);       

% Define cost functionx| and expected disturbances
% Q=[1300 0 0;0 1.8 0;0 0 1];
Q=[50 0 0;0 50 0;0 0 100000000];
R=[0.00001 0;0 0.0001];
W=ones(1,N)';  % expected demand (this is just an example)

[A B C]=model_system(vk,thetak,Ts);
[A,B,D,Q]=increase_matrixDUQ(A,B,D,Q);
Q(4,4)=5;
Q(5,5)=5;

[Gx,Gu,Gw]=constants_mpc(A,B,D,N);

% Build R_hat
R_hat = kron(eye(N),R);
% Build Q_hat
Q_hat=kron(eye(N),Q);

% Constraints
Ax=[1 0 0;0 1 0;-1 0 0 ;0 -1 0];
bx=100*[150; 150;150; 150];
Au=[1 0;0 1 ;-1 0;0 -1];
bu=[10; 4; 10; 4];
bux=bu;

Axaum=[Ax zeros(size(Ax,1),size(Au,2));zeros(size(Au,1),size(Ax,2)) Au];
bxaum=[bx;bu];
% Axaum=[Ax zeros(size(Ax,1),size(Au,2))];
% bxaum=[bx];
Ax=Axaum;
bx=bxaum;
bu=[8; 2; 8; 2];
% Transform into U constraints
Au_hat=kron(eye(N),Au);
bu_hat=kron(ones(N,1),bu);
Ax_hat=kron(eye(N),Ax);
bx_hat=kron(ones(N,1),bx);
%Delta U


% Aggregated U constraints
AU=[Ax_hat*Gu; Au_hat];
%bU=[bx_hat-Ax_hat*Gx*X0-Ax_hat*Gw*W;bu_hat];

% MPC into action
Xhist=X0;
Uhist=[];
VK=vk;
WK=wk;
THK=thetak;
Disturb= normrnd(0.5,1,Simlength+N,1); %Longer than simulation for prediction horizon
% Simulation loop
XR=[];
u=[0;0];
D=zeros(3,1);
path=createPath();
i=1;
delta=0.2;

figure();
plot(Xhist(1,:),Xhist(2,:),'b',path(1,:),path(2,:),'mo')
ylabel('y');
xlabel('X');
title('Trayectory');
grid on;
hold on;


for k=1:Simlength
    % expected disturbances (force that they are different)
%     W=-10*ones(N,1);%0*Disturb(k:k+N-1)+0*normrnd(0,0.2,N,1); 
    W=-10*Disturb(k:k+N-1)+0*normrnd(0,0.2,N,1);
    Wp=0*Disturb(k:k+N-1);
    % Update controller matrices for current state and disturbances (H and Au are constant)
    [A B C]=model_system(vk,thetak,Ts);
    [Xr,i]=createReferenceDU(path,i,X0,B,vk,Ts,N,delta);
    UMPC=MPC_DU(A,B,D,N,Wp,X0,Xr,Q_hat,R_hat,Au_hat,bu_hat,Ax_hat,bx_hat,u);
    % Apply only first component
    u=UMPC(1:size(B,2))+u;
%     X1=linearModel(A,B,D,u,Disturb(k),X0);
    X1=nonlinearModel(D,u,W(1),X0,thetak,wk,vk,Ts);
    vk=saturated(-bux(3),bux(1),vk+u(1));
    wk=saturated(-bux(4),bux(2),wk+u(2));
%     vk=(vk+u(1));
%     wk=(wk+u(2));
    thetak=X1(3);
    X0=X1;
    XR=[XR Xr(1:3)];
    VK=[VK vk];
    WK=[WK wk];
    THK=[THK thetak];
    Xhist=[Xhist X0];
    Uhist=[Uhist u];

%     plot(Xhist(1,:),Xhist(2,:),'b')
%     ylabel('y');
%     xlabel('X');
%     title('Trayectory');
%     grid on;
%     hold on;
%     pause(0.1);
%     [Gx,Gu,Gw]=constants_mpc(A,B,D,N);
%     Xpredict=Gx*X0+Gu*UMPC+Gw*W;
%     Xp=Xpredict(1:3:end);
%     Yp=Xpredict(2:3:end);
%     THp=Xpredict(3:3:end);
%     plot(Xp,Yp,'r');
%     v=[X0(1);X0(2)];
%     o= Xr(3);
%     m=2;
%     v=[v v+[m*cos(o);m*sin(o)]];
%     plot(v(1,:),v(2,:),'g');
    

end


%Simlength=size(Xhist,2);
t=0:Ts:Ts*(Simlength-1);
% figure();
% plot(Xhist);
% hold on;
% plot(Uhist);
% legend('X','U');
% xlim([0 Simlength]);

% state x and y 
figure();
plot(Xhist(1,:),Xhist(2,:),path(1,:),path(2,:),'o')
ylabel('y');
xlabel('X');
title('Trayectory');
grid on;

%u
figure();
subplot(4,1,1);
plot(t,Uhist(1,:));
ylabel('Velocity');
xlabel('Time(s)');
title('Delta Velocity');
grid on;

subplot(4,1,2);
plot(t,Uhist(2,:));
ylabel('Angular Velocity');
xlabel('Time(s)');
title('Delta Angular Velocity');
grid on;

% vk thetak
t2=0:Ts:Ts*(Simlength);
% figure();
subplot(4,1,3);
plot(t2,VK);
ylabel('Velocity');
xlabel('Time(s)');
title('Velocity');
grid on;

subplot(4,1,4);
plot(t2,WK);
%plot(THK);
ylabel('Angular Velocity');
xlabel('Time(s)');
title('Angular Velocity');
grid on;




t=0:Ts:Ts*(Simlength-1);
t2=0:Ts:Ts*(Simlength);
figure();
subplot(3,1,1);
plot(t,XR(1,:),t2,Xhist(1,:));
ylabel('Position');
xlabel('Time(s)');
title('Position on X axis');
grid on;

% vk thetak
t2=0:Ts:Ts*(Simlength);
% figure();
subplot(3,1,2);
plot(t,XR(2,:),t2,Xhist(2,:));
ylabel('Position');
xlabel('Time(s)');
title('Position on Y axis');
grid on;

subplot(3,1,3);
plot(t,XR(3,:),t2,THK);
%plot(THK);
ylabel('Theta');
xlabel('Time(s)');
title('Theta');
grid on;