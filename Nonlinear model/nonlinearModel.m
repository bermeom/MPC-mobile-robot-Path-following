function [X]=nonlinearModel(D,u,w,X0,thetak,wk,vk,Ts)
    X=X0+[(vk+u(1))*Ts*cos(X0(3));(vk+u(1))*Ts*sin(X0(3));(wk+u(2))*Ts]+D*w;
end