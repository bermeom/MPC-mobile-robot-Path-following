function [Xr,newi]=createReferenceDU(path,i,X0,B,vk,Ts,N,delta)
      n=1;
      r=vk/Ts;
      x=path(1,i);
      y=path(2,i);
      R=sqrt((x-X0(1))^2+(y-X0(2))^2);
      if(R<delta &&size(path,2)>i)
        i=i+1
      end
      j=i;
      Xr=[];
      while n<=N 
          x=path(1,j);
          y=path(2,j);
          R=sqrt((x-X0(1))^2+(y-X0(2))^2);
          t=R/vk;
         
          nn=ceil(t/Ts);
          if nn>N-n+1|| nn<0
             nn=N-n+1;
          end
          if size(path,2) ==j
               nn=N-n+1;
          end
%          nn=N-n+1;
          Xr=[Xr;kron(ones(nn,1),[x;y;(atan((y-X0(2))/(x-X0(1))));zeros(size(B,2),1)])];
          X0(1)=x;
          X0(2)=y;
          
          n=n+nn;
          j=j+1;
      end
     
      newi=i;
end