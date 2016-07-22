clc;
clear all;
close all;
% Initial state x(0)
X0=[0;0;(0*pi)/180];
vk=0;
Ts=0.01;
thetak=(0*pi)/180;
wk=0;
D=zeros(3,1);
N=10;
Xr=[1 1 0]';
Simlength=3*(1/Ts);

% Define cost function and expected disturbances
Q=[1 0 0;0 1 0;0 0 1];
R=[0.001 0;0 0.001];
W=ones(1,N)';  % expected demand (this is just an example)

[A B C]=model_system(vk,thetak,Ts);
[Gx,Gu,Gw]=constants_mpc(A,B,D,N);
% Build R_hat
R_hat = kron(eye(N),R);
% Build Q_hat
Q_hat=kron(eye(N),Q);


% Constraints
Ax=[1 0 0;0 1 0;-1 0 0 ;0 -1 0;0 0 1;0 0 -1];
bx=[150; 150;150; 150;1;1];
Au=[1 0;0 1;-1 0;0 -1];
bu=[150; 1; 0; 1];

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
for k=1:Simlength
    
    % expected disturbances (force that they are different)
    W=0*Disturb(k:k+N-1)+0*normrnd(0,0.2,N,1); 
    
    % Update controller matrices for current state and disturbances (H and Au are constant)
%     [A B C]=model_system(vk,thetak,Ts);
%     [Gx,Gu,Gw]=constants_mpc(A,B,D,N);
    % Build cost function
    H=Gu'*Q_hat*Gu+R_hat;
%     Xr(3)=(atan((Xr(2)-X0(2))/(Xr(1)-X0(1)))-Xr(3))/Ts;
    F=X0'*Gx'*Q_hat*Gu+W'*Gw'*Q_hat*Gu-kron(ones(N,1),Xr)'*Q_hat*Gu;
    % Aggregated U constraints
    AU=[Ax_hat*Gu; Au_hat];
    bU=[bx_hat-Ax_hat*Gx*X0-Ax_hat*Gw*W;bu_hat];
    UMPC=quadprog(H,F,AU,bU);
    XR=[XR Xr];
%     AU*UMPC2+bU
% Apply only first component
    if (size(UMPC)==[0 0])
        u=u2;
    else 
        u=UMPC(1:size(B,2));
        u2=UMPC(size(B,2)+1:(2)*size(B,2));
        UMPC2=UMPC;
        
    end
    
%     X1=A*X0+B*u+D*Disturb(k);
    X1=X0+[(vk+u(1))*Ts*cos(thetak);(vk+u(1))*Ts*sin(thetak);(wk+u(2))*Ts]+D*Disturb(k);
    
%     dx=(X1(1)-X0(1))/Ts;
%     dy=(X1(2)-X0(2))/Ts;
%     vk=sqrt(dx*dx+dy*dy);
    vk=vk+u(1);
    wk=wk+u(2);
    %thetak=atan(dy/dx);
    thetak=X1(3);
    VK=[VK vk];
    THK=[THK thetak];
    X0=X1;
    Xhist=[Xhist X0];
    Uhist=[Uhist u];
end
t=0:Ts:Ts*(Simlength-1);
figure();
plot(Xhist);
hold on;
plot(Uhist);
legend('X','U');
xlim([0 Simlength]);
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
plot(Xhist(2,:),Xhist(1,:))
ylabel('y');
xlabel('X');
title('Trayectory');
grid on;


%%
figure();
subplot(4,1,1);
plot(t,Uhist(1,:));
ylabel('Velocity');
xlabel('Time(s)');
title('Delta Velocity');
grid on;

subplot(3,1,1);
plot(t,XR(1,:));
ylabel('Velocity');
xlabel('Time(s)');
title('Velocity on X axis');
grid on;

% vk thetak
t2=0:Ts:Ts*(Simlength);
% figure();
subplot(3,1,2);
plot(t,XR(2,:));
ylabel('Velocity');
xlabel('Time(s)');
title('Velocity on Y axis');
grid on;

subplot(3,1,3);
plot(t,XR(3,:));
%plot(THK);
ylabel('Theta');
xlabel('Time(s)');
title('Theta');
grid on;
