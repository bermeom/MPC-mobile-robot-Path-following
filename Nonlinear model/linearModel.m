function [X]=linearModel(A,B,D,u,w,X0)
     X=A*X0+B*u+D*w;
end