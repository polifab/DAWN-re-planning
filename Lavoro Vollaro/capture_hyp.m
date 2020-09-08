function [traj2, delta_v2] = capture_hyp(goal_id, orbit, arrival_time,...
                             park_i, goal_coe, v2)
% ESCAPE_HYP computes the trajectory the spacecraft will follow
%   to escape desired planet's Sphere of Influence.
%
%   Options for the Sun are not contemplated since that would be
%   the general case of an interplanetary trajectory.
%
%   planet_id - identifier of the origin planet:
%                1 = Mercury
%                2 = Venus
%                3 = Earth
%                4 = Mars
%                5 = Jupiter
%                6 = Saturn
%                7 = Uranus
%                8 = Neptune
%                9 = Pluto
%               10 = Vesta
%               11 = Ceres
%
%   goal_id  - identifier of the destination planet:
%                1 = Mercury
%                2 = Venus
%                3 = Earth
%                4 = Mars
%                5 = Jupiter
%                6 = Saturn
%                7 = Uranus
%                8 = Neptune
%                9 = Pluto
%               10 = Vesta
%               11 = Ceres
%
%   orbit    - first two points of the interplanetary trajectory
%              computed via the patched conics method
%
%   arrival_time - array specifying time of arrivalarture with elements 
%                    (in this order):
%                     year         - range: 1901 - 2099
%                     month        - range: 1 - 12
%                     day          - range: 1 - 31
%                     hour         - range: 0 - 23
%                     minute       - range: 0 - 60
%                     second       - range: 0 - 60
%   park_r   - radius of the circular parking orbit around origin planet
%
%   goal_coe - classical orbital elements of the target interplanetary
%              orbit:
%                h    = angular momentum (km^2/s)
%                e    = eccentricity
%                RA   = right ascension of the ascending
%                       node (rad)
%                incl = inclination of the orbit (rad)
%                w    = argument of perigee (rad)
%                TA   = true anomaly (rad)
%                a    = semimajor axis (km)
%   v2 - velocity of the spacecraft at goal planet

    %% Argument validation
%     validateattributes(arrival_time,{'double'},{'size',[1 6]})
    %validateattributes(goal_coe,{'double'},{'size',[1 7]})

    %% Data
    global mu

    masses = 10^24 * [0.330
                      4.87
                      5.97
                      0.642
                      1898
                      568
                      86.8
                      102
                      0.0146
                      0.0002589
                      0.000947
                      1989100]; %[kg]

	radii = [2439.5
             6052 
             6378
             3396
             71492
             60268
             25559
             24764
             1185
             262.7
             476.2
             695508]; %[km] 

	distances = [57909227
                 108209475
                 149598262
                 227943824
                 778340821
                 1426666422
                 2870658186
                 4498396441
                 5906440628
                 491593189
                 423690250];%[km]
	
    aphelions = 10^6 * [57.91
                       108.21
                       149.60
                       227.92
                       778.57
                      1433.53
                      2872.46
                      4495.06
                      5906.38
                       353.35
                       414.087]; %[km]
    
    G    = 6.6742e-20; %[km^3/kg/s^2]
    
    %SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
    pl_SOI = (masses(goal_id)/masses(12))^(2/5)...
        * distances(goal_id); %[km]
    
    pl_mu = G * masses(goal_id); %[km^3/s^2]
    pl_radius = radii(goal_id); %[km]
    pl_a = aphelions(goal_id); %[km]
    goal_a = aphelions(goal_id); %[km]

    %% Needed variables
        [~, pl_r0, v_arrival, ~] =...
        planet_elements_and_sv(goal_id, arrival_time(1),arrival_time(2),...
                        arrival_time(3),arrival_time(4),arrival_time(5),arrival_time(6));
   
    V_arrival = sqrt(dot(v_arrival,v_arrival));       % module of v_arrival                   

    vinf= sqrt(dot(v2- v_arrival, v2 - v_arrival));          %sqrt(mu/pl_a)*(sqrt(2*goal_a/(pl_a+goal_a))-1);

    rp = pl_radius + park_i;
    e = 1 + rp*vinf^2/pl_mu;
    a = rp/(e-1);
    b = a*sqrt(e^2-1);
    Delta = sqrt(a*(1 - e^2)*pl_mu/vinf);    % raggio target per ottenere l'iperbole corretta
    half_delta = asin(1/e);                  % parametro iperbole
    
    v_hyp = sqrt(-2*pl_mu/rp + pl_mu/a);
    vc = sqrt(pl_mu/rp);                     % velocity of parking orbit
    delta_v2 = v_hyp - vc;                   % delta v from hyperbolic orbit to circular orbit

    beta = acos(1/e);
    

    h = Delta*vinf;
    RA = deg2rad(goal_coe(3));
    incl = deg2rad(goal_coe(4));
    w = deg2rad(goal_coe(5));

    n = sqrt(pl_mu/a^3);
    
    %% Trajectory computation
    rr = [];

    for t=0:60:24*3600
        M = n*t;
        F = kepler_H(e,M);
        cosf = (e-cosh(F))/(e*cosh(F)-1);
        f = acos(cosf);
        coe = [h, e, RA, incl, w, f];
        [r,~] = sv_from_coe(coe,pl_mu);
        rr = cat(1,rr,r);
    end

    %Angle of orientation of escape trajectory
    out_dir = orbit(2,1:3)-orbit(1,1:3);
    out_angle = deg2rad(atan2d_0_360(out_dir(2),out_dir(1)));

    t = 0:0.1:5;

    %Parametric hyperbola equations
    xh_l = -a*cosh(t);
    xh_r = a*cosh(t);
    yh = b*sinh(t);

    hyp = [];
    for i = 1:length(t)
        point = pl_r0' + Rotx(incl)*Rotx(park_i)*Rotz(out_angle)*Rotz(beta)*...
                    ([xh_r(i); -yh(i);0] + [-(a+rp);0;0]);
        hyp = cat(1,hyp,point');
        if norm(hyp(size(hyp,1),:)-hyp(1,:))>= pl_SOI
            break;
        end
    end
    
    %Hyperbola plot
    plot3(hyp(:,1),hyp(:,2),hyp(:,3),'mo-')
    grid on;
    
    traj = hyp;
end


                             