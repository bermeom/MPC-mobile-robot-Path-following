function [Aaum,Baum,Daum,Qaum]=increase_matrixDUQ(A,B,D,Q)
    Aaum=[A B;zeros(size(B,2),size(A,2)) eye(size(B,2))];
    Baum=[B;eye(size(B,2))];
    Daum=[D;zeros(size(B,2),size(D,2))];
    Qaum=[Q zeros(size(Q,2),size(B,2));zeros(size(B,2),size(Q,2)) eye(size(B,2))];
end 