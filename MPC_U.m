function [H,F,AU,bU]=MPC_U(A,B,D,N,W,X0,Xr,Q_hat,R_hat,Au_hat,bu_hat,Ax_hat,bx_hat)
    XR=referencePathGeneratorN(Xr,N);
    [Gx,Gu,Gw]=constants_mpc(A,B,D,N);
    % Build cost function
    H=Gu'*Q_hat*Gu+R_hat;
    F=X0'*Gx'*Q_hat*Gu+W'*Gw'*Q_hat*Gu-XR'*Q_hat*Gu;
    % Aggregated U constraints
    AU=[Ax_hat*Gu; Au_hat];
    bU=[bx_hat-Ax_hat*Gx*X0-Ax_hat*Gw*W;bu_hat];
end 