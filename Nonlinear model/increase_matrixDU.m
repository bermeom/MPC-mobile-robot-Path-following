function [Aaum,Baum,Daum]=increase_matrixDU(A,B,D)
    Aaum=[A B;zeros(size(B,2),size(A,2)) eye(size(B,2))];
    Baum=[B;eye(size(B,2))];
    Daum=[D;zeros(size(B,2),size(D,2))];
end 