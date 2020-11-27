
G = 6.67408e-20;   % universal gravitational constant [km^3/kg*s^2]  
m_sun = 1.989e+30; % solar mass kg
mu_sun = G*m_sun;  % standard gravitational parameter (sun) [km^3/s^2]
 
tf = 911 ;                % time of flight - days
conversion = 1.496e+8;    % km corresponding to 1 AU
r1_AU = [1.433816740218012 2.101176291368909 -2.376488907233350e-01]           %initial position in AU  05/09/12   
r2_AU = [2.330501934862029e-01 -2.854984767149892e+00 -1.322205442979184e-01]  %final position in AU    05/03/15   


r1 = r1_AU*conversion;     %initial position in km   05/09/12
r2 = r2_AU*conversion;     %final position in km     05/03/15

%% lambert_Rodyo
[V1, V2, extremal_distances, exitflag] = lambert_rodyo(r1, r2, tf, 0, mu_sun);

%V1 V2 = initial and final velocicties. [km/s]
%extremal distances =  minimum(1) and maximum(2) distance of the
%                             spacecraft to the central body.[km]

%% orbital elements

[a e i RAAN w f] = rv_2_orb(r1,V1,mu_sun)       % orbital parameters (refered to the spacecraft
[a1 e1 i1 RAAN1 w1 f1] = rv_2_orb(r2,V2,mu_sun) %                      trajectory from Vesta to Ceres)
% coe1 = coe_from_sv(r1,V1)                                        
% coe2 = coe_from_sv(r2,V2)



fprintf('\n  Eccentricity                                = %g',...
                                                           e)
fprintf('\n  Right ascension of the ascending node (rad) = %g',...
                                                           RAAN)
fprintf('\n  Inclination to the ecliptic (rad)           = %g',...
                                                           i)
fprintf('\n  Argument of perihelion (rad)                = %g',...
                                                           w)
fprintf('\n  True anomaly at departure (rad)             = %g',...
                                                           f)
fprintf('\n  True anomaly at arrival (rad)               = %g', ...
                                                           f1)
fprintf('\n  Semimajor axis (km)                         = %g\n',...
                                                           a)