clc;
clear all;
close all;
% Initial state x(0)
X0=[0;0;(0*pi)/180];
%X0 =[2.3917 ;-1.3973;-1.2586];
vk=0;
Ts=0.005;
%thetak=(120*pi)/180;
thetak=-1.2586;
wk=0;
D=zeros(3,1);
N=50;
Xr=[-2 10 atan(1)]';
% Xr(3)=(atan((-50-X0(2))/(50-X0(1))));
% Xr(3)=(atan((Xr(2)-X0(2))/(Xr(1)-X0(1))));
s=3;
Simlength=(s/Ts);       

% Define cost functionx| and expected disturbances
Q=[0.00001 0 0;0 1 0;0 0 100000000];
R=[0.01 0;0 100];
W=ones(1,N)';  % expected demand (this is just an example)

[A B C]=model_system(vk,thetak,Ts);
[Gx,Gu,Gw]=constants_mpc(A,B,D,N);
% Build R_hat
R_hat = kron(eye(N),R);
% Build Q_hat
Q_hat=kron(eye(N),Q);

% Constraints
Ax=[1 0 0;0 1 0;-1 0 0 ;0 -1 0];
bx=100*[150; 150;150; 150];
Au=[1 0;0 1 ;-1 0;0 -1];
bu=[2; 0.01; 0; 0.01];

% Transform into U constraints
Au_hat=kron(eye(N),Au);
bu_hat=kron(ones(N,1),bu);
Ax_hat=kron(eye(N),Ax);
bx_hat=kron(ones(N,1),bx);

% Aggregated U constraints
AU=[Ax_hat*Gu; Au_hat];
bU=[bx_hat-Ax_hat*Gx*X0-Ax_hat*Gw*W;bu_hat];

% MPC into action
Xhist=X0;
Uhist=[];
VK=vk;
THK=thetak;
Disturb= normrnd(0.5,1,Simlength+N,1); %Longer than simulation for prediction horizon
% Simulation loop
XR=[];
% figure();
% plot(Xhist(1,:),Xhist(2,:),'b')
% ylabel('y');
% xlabel('X');
% title('Trayectory');
% grid on;
% hold on;
for k=1:Simlength
    % expected disturbances (force that they are different)
    W=0*Disturb(k:k+N-1)+0*normrnd(0,0.2,N,1); 
    % Update controller matrices for current state and disturbances (H and Au are constant)
    [A B C]=model_system(vk,thetak,Ts);
    Xr(3)=(atan((Xr(2)-X0(2))/(Xr(1)-X0(1))));
    %Xr(3)=(atan((-50-X0(2))/(50-X0(1))));
    [H,F,AU,bU]=MPC_U(A,B,D,N,W,X0,Xr,Q_hat,R_hat,Au_hat,bu_hat,Ax_hat,bx_hat);
    UMPC=quadprog(2*H,2*F,AU,bU);
    XR=[XR Xr];
    % Apply only first component
    u=UMPC(1:size(B,2));
%     X1=linearModel(A,B,D,u,Disturb(k),X0);
    X1=nonlinearModel(D,u,Disturb(k),X0,thetak,wk,vk,Ts);

    vk=saturated(0,20,vk+u(1));
    wk=saturated(-1,1,wk+u(2));
    thetak=X1(3);
    X0=X1;
    VK=[VK vk];
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
%     Xpredict=Gx*X0+Gu*UMPC+Gw*W;
%     Xp=Xpredict(1:3:end);
%     Yp=Xpredict(2:3:end);
%     THp=Xpredict(3:3:end);
%     plot(Xp,Yp,'R')

end

%%

%Simlength=size(Xhist,2);
t=0:Ts:Ts*(Simlength-1);
% figure();
% plot(Xhist);
% hold on;
% plot(Uhist);
% legend('X','U');
% xlim([0 Simlength]);
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
plot(t2,THK,t2,Xhist(3,:));
%plot(THK);
ylabel('Theta');
xlabel('Time(s)');
title('Theta');
grid on;

% state x and y 
figure();
plot(Xhist(1,:),Xhist(2,:))
ylabel('y');
xlabel('X');
title('Trayectory');
grid on;


%%
t=0:Ts:Ts*(Simlength-1);
t2=0:Ts:Ts*(Simlength);
figure();
subplot(3,1,1);
plot(t,XR(1,:),t2,Xhist(1,:));
ylabel('Velocity');
xlabel('Time(s)');
title('Velocity on X axis');
grid on;

% vk thetak
t2=0:Ts:Ts*(Simlength);
% figure();
subplot(3,1,2);
plot(t,XR(2,:),t2,Xhist(2,:));
ylabel('Velocity');
xlabel('Time(s)');
title('Velocity on Y axis');
grid on;

subplot(3,1,3);
plot(t,XR(3,:),t2,THK);
%plot(THK);
ylabel('Theta');
xlabel('Time(s)');
title('Theta');
grid on;