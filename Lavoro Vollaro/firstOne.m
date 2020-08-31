%%                      WARNING
%       Collection of all code used so far:
%       some of it is purposeless, some is garbage!



%% Data

%Some of these can be derived from planet_elements_and_sv
Earth_ecc = 0.0167;
Earth_major = 149.6*10^6; %[km] <-> 1AU
Earth_minor = Earth_major*sqrt(1-Earth_ecc^2); %[km]
Earth_peri = 147.09*10^6; %[km]
Earth_radius = 6378; %[km]
Earth_c = sqrt(Earth_major^2-Earth_minor^2); %[km]

Mars_ecc = 0.0935;
Mars_major = 227.92*10^6; %[km]
Mars_minor = Mars_major*sqrt(1-Mars_ecc^2); %[km]
Mars_peri = 206.62*10^6; %[km]
Mars_radius = 3390; %[km]
Mars_c = sqrt(Mars_major^2-Mars_minor^2); %[km]

Ceres_to_Sun = 446145795; %[km]
Ceres_mass = 8.958*10^20; %[kg]

Vesta_to_Sun = 353*10^6; %[km]
Vesta_mass = 2.589*10^20; %[kg]

Sun_mass = 1.989*10^30; %[kg]
Sun_radius = 696.340; %[km]

%SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
Earth_SOI = 9.24*10^5; %[km]
Mars_SOI = 5.74*10^5; %[km]
Ceres_SOI = (Ceres_mass/Sun_mass)^(2/5)*Ceres_to_Sun; %[km]
Vesta_SOI = (Vesta_mass/Sun_mass)^(2/5)*Vesta_to_Sun; %[km]

parkingOrbit = Earth_radius+200; %[km]

Earth_mu = 398600; %[km^3/s^2]
Mars_mu = 42828;%[km^3/s^2]

%Sun gravitational parameter
global mu
mu = 1.327*(10^11); %[km^3/s^2]

%Park orbit
Park_r0 = [parkingOrbit,0,0]; %TO BE CHANGED AS NEEDED
%for circular orbits: 
%   v = sqrt(mu/r)
Park_v0 = sqrt(Earth_mu/Park_r0(1)); %km/s

%% Initial conditions
%Parking orbit right ascension and declination
[Park_ra, Park_dec] = ra_and_dec_from_r(Park_r0);

%Earth elements and sv in the eliocentric system
% with respect to the sun
[Earth_coe0, Earth_r0, Earth_v0, Earth_julianday] = planet_elements_and_sv(3,2007,9,27,0,0,0);

%Mars elements and sv in the eliocentric system
% with respect to the sun
[Mars_coe0, Mars_r0, Mars_v0, Mars_julianday] = planet_elements_and_sv(4,2007,9,27,0,0,0);

%Mars right ascension and declination
[Mars_ra,Mars_dec] = ra_and_dec_from_r(Mars_r0);

%Angle between Earth and Mars
EM_cos = [Earth_r0(1:2),0]*[Mars_r0(1:2),0]'/(norm([Earth_r0(1:2),0])*norm([Mars_r0(1:2),0]));
Earth_Mars = rad2deg(acos(EM_cos));

%% Final conditions
[Earth_coef, Earth_rf, Earth_vf, Earth_juliandayf] = planet_elements_and_sv(3,2009,2,17,0,0,0);
[Mars_coef, Mars_rf, Mars_vf, Mars_juliandayf] = planet_elements_and_sv(4,2009,2,17,0,0,0);

%% Foci
Earth_focus = Rotz(Earth_coe0(3))*Rotx(Earth_coe0(4))*[0;sqrt(Earth_coe0(7)^2-(Earth_coe0(7)*sqrt(1-Earth_coe0(2)^2))^2);0];
Mars_focus = Rotz(Mars_coe0(3))*Rotx(Mars_coe0(4))*[0;sqrt(Mars_coe0(7)^2-(Mars_coe0(7)*sqrt(1-Mars_coe0(2)^2))^2);0];

%% Circular orbit
%Period
Park_T = 2*pi*parkingOrbit^(3/2)/sqrt(Earth_mu);

%Time between two positions:
%   t = phi/(2*pi)*Park_T;
% t_toEarthPerigee = deg2rad(angle_toPerigee)/(2*pi)*Park_T; %[s]

%% For testing purposes
figure(1)
grid
hold on
Earth_park(Earth_r0)
orbit = orbitAttempt;

%% Choosing point of departure
n = 10^10;
day = 0;

for g = 1:687
    %New planets' positions
    [E_r, E_v] = rv_from_r0v0(Earth_r0, Earth_v0, g*60*60*24);
    [M_r, M_v] = rv_from_r0v0(Mars_r0, Mars_v0, g*60*60*24);
    EM_c = [E_r(1:2),0]*[M_r(1:2),0]'/(norm([E_r(1:2),0])*norm([M_r(1:2),0]));
    
    %To find the perihelium
    if(norm(E_r)<n)
        day = g;
        n = norm(E_r);
    end
    
    %Angle between planets (might be wrong atm)
    E_coe = coe_from_sv(E_r, E_v, Earth_mu);
    M_coe = coe_from_sv(M_r, M_v, Mars_mu);
    diff = E_r-M_r;
    if diff(1)==0
        EM_angle = -sign(M_r(1))*sign(diff(1))*rad2deg(acos(EM_c));
    else
        EM_angle = sign(diff(1))*rad2deg(acos(EM_c));
    end

    %to visualize the position
%     pp = Rotx(Earth_coe0(4))'*Rotz(Earth_coe0(3))'*E_r';
%     Rotx(Mars_coe0(4))'*Rotz(Mars_coe0(3))'*M_r';
%     plot3(E_r(1),E_r(2),E_r(3),'b+')
%     plot3(M_r(1),M_r(2),M_r(3),'ro')
%     plot3(pp(1),pp(2),pp(3),'b.')

    %the angle is usually around 44°,
    %minus safety margin for further operations

    %to find the moment in which the angle is as desired
%     if(EM_angle >= 0) && (EM_angle < 45)%?
%         fprintf('The angle between Earth and Mars is %4.2f after %d days \n',EM_angle,g)
%         break;
%     end
end

%% After some time
%'day' indicates Earth's perihelium

% day = day + 365;

%Earth (r,v) from (r0,v0)
[Earth_r, Earth_v] = rv_from_r0v0(Earth_r0, Earth_v0, day*60*60*24);

%Mars (r,v) from (r0,v0)
[Mars_r, Mars_v] = rv_from_r0v0(Mars_r0, Mars_v0, day*60*60*24);

%Mars right ascension and declination
[Mars_ra2,Mars_dec2] = ra_and_dec_from_r(Mars_r);

%Earth classical orbital elements from state vector
Earth_coe = coe_from_sv(Earth_r, Earth_v, Earth_mu);

%Mars classical orbital elements from state vector
Mars_coe = coe_from_sv(Mars_r, Mars_v, Mars_mu);

%% Orbit plots
alpha = 0:pi/200:2*pi;
phi = -pi/2:pi/100:pi/2;

plot_orbit(Earth_coe0(3),Earth_coe0(4),Earth_major,Earth_minor,Earth_focus)
xlabel('x')
ylabel('y')
zlabel('z')

plot_orbit(Mars_coe0(3),Mars_coe0(4),Mars_major,Mars_minor,Mars_focus)

plot3(Earth_r0(1),Earth_r0(2),Earth_r0(3),'bo')
plot3(Earth_rf(1),Earth_rf(2),Earth_rf(3),'bx')
plot3(Mars_r0(1),Mars_r0(2),Mars_r0(3),'ro')
plot3(Mars_rf(1),Mars_rf(2),Mars_rf(3),'rx')

% plot3([Earth_focus(1),Earth_r(1)],[Earth_focus(2),Earth_r(2)],[Earth_focus(3),Earth_r(3)],'b-+')
% plot3([Mars_focus(1),Mars_r(1)],[Mars_focus(2),Mars_r(2)],[Mars_focus(3),Mars_r(3)],'r-+')

%% Lambert attempt
part = datevec(datenum(2007,12,18));
arr = datevec(datenum(2009,2,17));
diff = etime(arr,part);
E = planet_elements_and_sv(3,2007,12,18,0,0,0);
M = planet_elements_and_sv(4,2009,2,27,0,0,0);