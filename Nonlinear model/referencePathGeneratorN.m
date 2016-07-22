function [XR]=ReferencePathGenerator(Xr,N)
    XR=kron(ones(N,1),Xr);
end 