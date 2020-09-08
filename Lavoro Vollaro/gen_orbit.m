function [dep_r,dep_v,arr_r,arr_v,flight,orb_oe] = ...
                        gen_orbit(dep_id,arr_id, dep_time, arr_time)
% GEN_ORBIT(dep_id, arr_id, dep_time, arr_time) generates info about
%   the interplanetary orbit of a 1000kg spacecraft departing from 
%   DEP_ID at DEP_TIME and arriving to ARR_ID at ARR_TIME.
%   [dep_r,dep_v,arr_r,arr_v,flight,orb_oe] = GEN_ORBIT(...)
%   returns the starting and final position DEP_R, ARR_R
%   of the spacecraft (corresponding to the main bodies position
%   due to the use of the patched conics method), along with its
%   starting and final velocities DEP_V, ARR_V, the time of flight
%   FLIGHT and the orbital elements ORB_OE of the trajectory.
%   
%    dep_id,arr_id - departure/arrival body identifier:
%                   1 = Mercury
%                   2 = Venus
%                   3 = Earth
%                   4 = Mars
%                   5 = Jupiter
%                   6 = Saturn
%                   7 = Uranus
%                   8 = Neptune
%                   9 = Pluto
%                  10 = Vesta
%                  11 = Ceres
%                  12 = Sun
%
%     dep_time     - array specifying time of departure with elements 
%                    (in this order):
%                     year         - range: 1901 - 2099
%                     month        - range: 1 - 12
%                     day          - range: 1 - 31
%                     hour         - range: 0 - 23
%                     minute       - range: 0 - 60
%                     second       - range: 0 - 60 
% 
%     arr_time     - array specifying time of arrival with elements 
%                    (in this order):
%                     year         - range: 1901 - 2099
%                     month        - range: 1 - 12
%                     day          - range: 1 - 31
%                     hour         - range: 0 - 23
%                     minute       - range: 0 - 60
%                     second       - range: 0 - 60 
%
%     orb_oe - orbital elements [h e RA incl w TA a]
%                     where
%                    h    = angular momentum (km^2/s)
%                    e    = eccentricity
%                    RA   = right ascension of the ascending
%                           node (rad)
%                    incl = inclination of the orbit (rad)
%                    w    = argument of perigee (rad)
%                    TA   = true anomaly (rad)
%                    a    = semimajor axis (km)
% 
%     flight          - time of flight from planet 1 to planet 2 (days)
% 
% 
%     vinf1, vinf2 - hyperbolic excess velocities at departure
%                  and arrival (km/s)
%     dep_r        - position of the main body at departure
%     dep_v        - velocity of the spacecraft at departure from
%                    origin SOI
%     arr_r        - position of the main body at arrival
%     arr_v        - velocity of the spacecraft at arrival in the 
%                    destination SOI
% 
%     User M-functions required: interplanetary, coe_from_sv,
%                              month_planet_names

    %% Argument validation
    validateattributes(dep_time,{'double'},{'size',[1 6]})
    validateattributes(arr_time,{'double'},{'size',[1 6]})

    %% Data
    %Sun mu
    global mu
    deg = pi/180;

    % Departure
    planet_id = dep_id;
    year      = dep_time(1);
    month     = dep_time(2);
    day       = dep_time(3);
    hour      = dep_time(4);
    minute    = dep_time(5);
    second    = dep_time(6);
    depart = [planet_id  year  month  day  hour  minute  second];

    % Arrival
    planet_id = arr_id;
    year      = arr_time(1);
    month     = arr_time(2);
    day       = arr_time(3);
    hour      = arr_time(4);
    minute    = arr_time(5);
    second    = arr_time(6);
    arrive = [planet_id  year  month  day  hour  minute  second];

    %% Computations
    [planet1, planet2, trajectory] = interplanetary(depart, arrive);

    %Rp1, Vp1: state vector of planet 1 at departure (km, km/s)
    %R1, V1: heliocentric state vector of the spacecraft 
    %at departure (km, km/s)
    R1  = planet1(1,1:3);
    Vp1 = planet1(1,4:6);
    jd1 = planet1(1,7);

    %Rp2, Vp2: state vector of planet 2 at arrival (km, km/s)
    %R2, V2: heliocentric state vector of the spacecraft at
    %arrival (km, km/s)
    R2  = planet2(1,1:3);
    Vp2 = planet2(1,4:6);
    jd2 = planet2(1,7);

    V1  = trajectory(1,1:3);
    V2  = trajectory(1,4:6);

    tof = jd2 - jd1;

    % Use Algorithm 4.2 to find the orbital elements of the
    % spacecraft trajectory based on [Rp1, V1]...
    coe  = coe_from_sv(R1, V1, mu);
    %   ... and [R2, V2]
    coe2 = coe_from_sv(R2, V2, mu);

    % Equations 8.94 and 8.95:
    vinf1 = V1 - Vp1;
    vinf2 = V2 - Vp2;

    % Echo the input data and output the solution to
    % the command window:
    [mm, pp] = month_planet_names(depart(3),depart(1));
    
    %% Output info
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('-----------------------------------------------------')
    fprintf('\n\n Departure:\n');
    fprintf('\n   Planet: %s', pp);%planet_name(depart(1)))
    fprintf('\n   Year  : %g', depart(2))
    fprintf('\n   Month : %s', mm);%month_name(depart(3)))
    fprintf('\n   Day   : %g', depart(4))
    fprintf('\n   Hour  : %g', depart(5))
    fprintf('\n   Minute: %g', depart(6))
    fprintf('\n   Second: %g', depart(7))
    fprintf('\n   Planet position vector (km)    = [%g  %g  %g]', ...
                                                   R1(1),R1(2), R1(3))

    fprintf('\n   Magnitude                      = %g\n', norm(R1))

    fprintf('\n   Planet velocity (km/s)         = [%g  %g  %g]', ...
                                     Vp1(1), Vp1(2), Vp1(3))

    fprintf('\n   Magnitude                      = %g\n', norm(Vp1))

    fprintf('\n   Spacecraft velocity (km/s)     = [%g  %g  %g]', ...
                                                   V1(1), V1(2), V1(3))

    fprintf('\n   Magnitude                      = %g\n', norm(V1))

    fprintf('\n   v-infinity at departure (km/s) = [%g  %g  %g]', ...
                                           vinf1(1), vinf1(2), vinf1(3))

    fprintf('\n   Magnitude                      = %g\n', norm(vinf1))
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n\n  ~~~~~~~~~~~~\n')
    fprintf('\n\n Time of flight = %g days\n', tof)
    fprintf('\n\n  ~~~~~~~~~~~~\n')
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [mm,pp] = month_planet_names(arrive(3),arrive(1));
    fprintf('\n\n Arrival:\n');
    fprintf('\n   Planet: %s', pp);
    fprintf('\n   Year  : %g', arrive(2))
    fprintf('\n   Month : %s', mm);
    fprintf('\n   Day   : %g', arrive(4))
    fprintf('\n   Hour  : %g', arrive(5))
    fprintf('\n   Minute: %g', arrive(6))
    fprintf('\n   Second: %g', arrive(7))
    fprintf('\n   Planet position vector (km)   = [%g  %g  %g]', ...
                                                  R2(1), R2(2), R2(3))

    fprintf('\n   Magnitude                     = %g\n', norm(R1))

    fprintf('\n   Planet velocity (km/s)        = [%g  %g  %g]', ...
                                      Vp2(1), Vp2(2), Vp2(3))

    fprintf('\n   Magnitude                     = %g\n', norm(Vp2))

    fprintf('\n   Spacecraft Velocity (km/s)    = [%g  %g  %g]', ...
                                                  V2(1), V2(2), V2(3))

    fprintf('\n   Magnitude                     = %g\n', norm(V2))

    fprintf('\n   v-infinity at arrival (km/s)  = [%g  %g  %g]', ...
                                         vinf2(1), vinf2(2), vinf2(3))

    fprintf('\n   Magnitude                     = %g', norm(vinf2))

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n\n  ~~~~~~~~~~~~\n')
    fprintf('\n\n\n Orbital elements of flight trajectory:\n')

    fprintf('\n  Angular momentum (km^2/s)                   = %g',...
                                                               coe(1))
    fprintf('\n  Eccentricity                                = %g',...
                                                               coe(2))
    fprintf('\n  Right ascension of the ascending node (deg) = %g',...
                                                           coe(3)/deg)
    fprintf('\n  Inclination to the ecliptic (deg)           = %g',...
                                                           coe(4)/deg)
    fprintf('\n  Argument of perihelion (deg)                = %g',...
                                                           coe(5)/deg)
    fprintf('\n  True anomaly at departure (deg)             = %g',...
                                                           coe(6)/deg)
    fprintf('\n  True anomaly at arrival (deg)               = %g\n', ...
                                                          coe2(6)/deg)
    fprintf('\n  Semimajor axis (km)                         = %g',...
                                                               coe(7))
	%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % If the orbit is an ellipse, output the period:
    if coe(2) < 1
        fprintf('\n  Period (days)                               = %g', ...
                                          2*pi/sqrt(mu)*coe(7)^1.5/24/3600)
    end
    fprintf('\n-----------------------------------------------------\n')
    
    %% Output arguments
    dep_r = R1;
    dep_v = V1;
    arr_r = R2;
    arr_v = V2;
    flight = tof;
    %         h     , e     , RA    , incl  , w     , TA    , a
    orb_oe = [coe(1), coe(2), coe(3), coe(4), coe(5), coe(6), coe(7)];
end