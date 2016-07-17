clc;
clear all;
close all;
% Initial state x(0)
X0=[1;0;0.01];
vk=1;
Ts=0.2;
thetak=0.5;
D=zeros(3,1);
N=3;
Xr=[50 50 0]';

% Define cost function and expected disturbances
Q=eye(3);
R=eye(2);
W=ones(1,N)';  % expected demand (this is just an example)

[A B C]=model_system(vk,thetak,Ts);
[Gx,Gu,Gw]=constants_mpc(A,B,D,N);
% Build R_hat
R_hat = kron(eye(N),R);
% Build Q_hat
Q_hat=kron(eye(N),Q);


% Constraints
Ax=[1 0 0;0 1 0;-1 0 0 ;0 -1 0];
bx=[50; 50;50; 50];
Au=[1 0;0 1;-1 0;0 -1];
bu=[150; 1; 25; 1];

% Transform into U constraints
Au_hat=kron(eye(N),Au);
bu_hat=kron(ones(N,1),bu);
Ax_hat=kron(eye(N),Ax);
bx_hat=kron(ones(N,1),bx);

% Aggregated U constraints
AU=[Ax_hat*Gu; Au_hat];
bU=[bx_hat-Ax_hat*Gx*X0-Ax_hat*Gw*W;bu_hat];

% MPC into action
Simlength=100;
Xhist=X0;
Uhist=[];
Disturb= normrnd(0.5,1,Simlength+N,1); %Longer than simulation for prediction horizon
% Simulation loop
for k=1:Simlength
    
    % expected disturbances (force that they are different)
    W=Disturb(k:k+N-1)+0*normrnd(0,0.2,N,1); 
    
    % Update controller matrices for current state and disturbances (H and Au are constant)
    [A B C]=model_system(vk,thetak,Ts);
    [Gx,Gu,Gw]=constants_mpc(A,B,D,N);
    % Build cost function
    H=Gu'*Q_hat*Gu+R_hat;
    
    F=X0'*Gx'*Q_hat*Gu+W'*Gw'*Q_hat*Gu-kron(ones(N,1),Xr)'*Q_hat*Gu;
    %bU=[bx_hat-Ax_hat*Gx*x-Ax_hat*Gw*W;bu_hat];
    
    UMPC=quadprog(H,F,AU,bU);
    %UMPC=quadprog(H,F);
    
    % Apply only first component
    u=UMPC(1:size(B,2));
    X1=A*X0+B*u+D*Disturb(k);
    dx=(X1(1)-X0(1))/Ts;
    dy=(X1(2)-X0(2))/Ts;
    %vk=sqrt(dx*dx+dy*dy);
    %thetak=arctan(dy/dx);
    X0=X1;
    Xhist=[Xhist X0];
    Uhist=[Uhist u];
end

plot(Xhist);
hold on;
plot(Uhist);
legend('X','U');
xlim([0 Simlength]);
figure();
plot(Xhist(1,:),Xhist(2,:))



