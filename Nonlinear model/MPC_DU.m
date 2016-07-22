function [UMPC]=MPC_DU(A,B,D,N,W,X0,Xr,Q_hat,R_hat,Au_hat,bu_hat,Ax_hat,bx_hat,um1)
    [Aaum,Baum,Daum]=increase_matrixDU(A,B,D);
    X0=[X0;um1];
    [H,F,AU,bU]=MPC_U2(Aaum,Baum,Daum,N,W,X0,Xr,Q_hat,R_hat,Au_hat,bu_hat,Ax_hat,bx_hat);
    UMPC=quadprog(2*H,F,AU,bU);
end 