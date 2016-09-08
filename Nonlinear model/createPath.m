function [path]=createPath()

    x = 0:20;
%     y = 30*sin(x/5);
    y = 30*cos(x/5);
%     y = 30*cos(x/5+10)-60*sin(x/5);
    xx = 1:1:20;
    yy = spline(x,y,xx);
    
    path=[xx;yy];
%     path=[xx(1);yy(1)];
%     path=[path [xx(2);yy(2)]];

%     figure();
%     plot(xx,yy,'o')
%     ylabel('y');
%     xlabel('X');
%     title('Trayectory');
%     grid on;
%    
end