function plot_orbit(ra,incl,maj,min,focus)
    %argument of cosine/sine
    alpha = 0:pi/100:2*pi;
    
    %rotation matrix
    rx = Rotx(incl);
    rz = Rotz(ra);
    rr = Rotz(deg2rad(48.7988));
    
    %standard ellipse equations
    x = maj*cos(alpha);
    y = min*sin(alpha);
    z = zeros(length(x));
    
    i = 1:1:length(x);
    
    %Rotated ellipse
    points  = rz*rx*rr*[x(i);y(i);z(i)];
    
    plot3(focus(1)+points(1,1:size(points,2)),focus(2)+points(2,1:size(points,2)),focus(3)+points(3,1:size(points,2)),'--')
    
    %tested alternatives
%     plot3(-0.05*10^8+points(1,1:size(points,2)),0.17*10^8+2.4983*10^6+points(2,1:size(points,2)),points(3,1:size(points,2))) -> Mars
%     plot3(points(1,1:size(points,2)),-3.0*10^6+points(2,1:size(points,2)),points(3,1:size(points,2)))
%     plot3(points(1,1:size(points,2)),points(2,1:size(points,2)),points(3,1:size(points,2)))
end