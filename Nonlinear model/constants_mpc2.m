function [Gx,Gu,Gw]=constants_mpc2(A,B,D,N)
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



