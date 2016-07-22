function [A,B,C] =model_system(vk,thetak,Ts)
    A=[1 0 -vk*Ts*sin(thetak);0 1 vk*Ts*cos(thetak); 0 0 1];
    B=[Ts*cos(thetak) -0.5*Ts*Ts*vk*sin(thetak);Ts*sin(thetak) 0.5*Ts*Ts*vk*cos(thetak);0 Ts];
    C=A(1:2,:);
end