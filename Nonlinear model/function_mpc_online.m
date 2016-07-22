function y=function_mpc_online(u)
global X0 vk thetak N Ts Q R D W
[A,B,C]=model_system(vk,thetak,Ts);

% Define prediction horizon and other useful variables
nx=size(A,2); % number of states
nu=size(B,2); % number of inputs
nw=size(D,2); % number of disturbances

% Future states can be written as X=Gx*x(0)+Gu*U+Gw*W
% Build Gx
% Build Gu
Gu=zeros(N*nx,N*nu);
% Build Gw
Gw=zeros(N*nx,N*nw);
for i=1:N
   Gx((i-1)*nx+1:i*nx,:)=A^i;
   for j=1:i
      Gu((i-1)*nx+1:(i)*nx,(j-1)*nu+1:j*nu)=A^(i-j)*B;
      Gw((i-1)*nx+1:(i)*nx,(j-1)*nw+1:j*nw)=A^(i-j)*D;
   end
end

% Build R_hat
R_hat = kron(eye(N),R);

% Build Q_hat
Q_hat=kron(eye(N),Q);


% Build cost function
H=Gu'*Q_hat*Gu+R_hat;
F=X0'*Gx'*Q_hat*Gu+W'*Gw'*Q_hat'*Gu;

%UMPC=quadprog(H,F,AU,bU);
UMPC=quadprog(H,F);
    
% Apply only first component
u=UMPC(1);
x=A*X0+B*u+D*Disturb(k);
Xhist=[Xhist x];
Uhist=[Uhist u];



