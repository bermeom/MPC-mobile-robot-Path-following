function [j]=binarySearch(vector,data,delta,i,l,r)
        j=-1;
        while (j==-1&&(vector(l)<=data&&vector(r)>=data))
               m=(l+r)/2;
                if (vector(m)-delta<data)
        end
        

end