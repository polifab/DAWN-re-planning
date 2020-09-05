% clc
% clear all
% close all

global mu
% mu = 1.327*(10^11); %[km^3/s^2]

[Earth_coe0, Earth_r0, Earth_v0, Earth_julianday] =...
    planet_elements_and_sv(3,2007,9,27,0,0,0);

Earth_SOI = 9.24*10^5; %[km]
Earth_mu = 398600; %[km^3/s^2]
Earth_radius = 6378; %[km]
Earth_a = 149.6*10^6; %[km]
Mars_a = 227.92*10^6; %[km]

vinf = sqrt(mu/Earth_a)*(sqrt(2*Mars_a/(Earth_a+Mars_a))-1);

rp = Earth_radius+200;
e = 1+rp*vinf^2/Earth_mu;
a = rp/(e-1);

vc = sqrt(Earth_mu/rp);
vp = sqrt(vinf^2+2*Earth_mu/rp);

beta = acos(1/e);

h = rp*vp;
RA = deg2rad(3.41309);
incl = deg2rad(1.84517);
w = deg2rad(45.1022);

n = sqrt(Earth_mu/a^3);
% figure2()
% hold on
% grid

% park_orbit(3,Earth_r0,200)
rr = [];

for t=0:60:24*3600
    M = n*t;
    F = kepler_H(e,M);
    cosf = (e-cosh(F))/(e*cosh(F)-1);
    f = acos(cosf);
%     fprintf("TA: %4.2f \n",rad2deg(f))
    coe = [h, e, RA, incl, w, f];
    [r,v] = sv_from_coe(coe,Earth_mu);
    rr = cat(1,rr,r);
%     plot3(Earth_r0(1)+r(1),Earth_r0(2)+r(2),Earth_r0(3)+r(3),'bx')
end
% [xx,yy,zz] = sphere(10);
% surface(Earth_r0(1)+Earth_SOI*xx, Earth_r0(2)+Earth_SOI*yy,...
%         Earth_r0(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor','b')

% xlabel('x')
% ylabel('y')
% zlabel('z')
% xlim([ 1.486*10^8, 1.508*10^8])
% ylim([8*10^6, 10*10^6])
% zlim([-10*10^5, 10*10^5])
% view(-10,45)
plot3(Earth_r0(1)+rr(:,1),Earth_r0(2)+rr(:,2),Earth_r0(3)+rr(:,3),'bo-')

%Angle of orientation of Earth velocity
vel_dir = deg2rad(atan2d_0_360(Earth_v0(2),Earth_v0(1)));

%To visualize Earth velocity direction
vect = [Earth_radius+200;0;0];
v_aligned = Earth_r0' + Rotz(vel_dir)*vect;
plot3([Earth_r0(1),v_aligned(1)],[Earth_r0(2),v_aligned(2)],...
    [Earth_r0(3),v_aligned(3)],'ro-')
vect = [-(Earth_radius+200);0;0];
v_aligned = Earth_r0' + Rotz(vel_dir)*vect;
plot3([Earth_r0(1),v_aligned(1)],[Earth_r0(2),v_aligned(2)],...
    [Earth_r0(3),v_aligned(3)],'ro-')

c_hyp = Earth_r0' + Rotz(vel_dir)*Rotz(beta)*[-(Earth_radius+200);0;0];
plot3([Earth_r0(1),c_hyp(1)],[Earth_r0(2),c_hyp(2)],...
    [Earth_r0(3),c_hyp(3)],'go-')
