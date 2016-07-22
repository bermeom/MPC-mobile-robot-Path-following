function [Xr,newi]=createReference(path,i,X0,vk,Ts,N,delta)
      n=1;
      r=vk/Ts;
      x=path(1,i);
      y=path(2,i);
      R=sqrt((x-X0(1))^2+(y-X0(2))^2);
      if(R<delta)
        i=i+1;
      end
      j=i
      Xr=[];
      XR=kron(ones(N,1),Xr)
      while n<N 
          n
          x=path(1,j);
          y=path(2,j);
          R=sqrt((x-X0(1))^2+(y-X0(2))^2);
          t=R/vk;
          ceil(t/Ts)
          N-n+1
          nn=min(ceil(t/Ts),N-n+1)
          if size(path,2) ==j
            nn=N-n+1;
          end
          Xr=[Xr;kron(ones(nn,1),[x;y;(atan((y-X0(2))/(x-X0(1))))])]
          n=n+nn
          j=j+1;
      end
      newi=i;
end