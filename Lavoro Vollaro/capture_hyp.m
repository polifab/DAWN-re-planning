function [traj, delta_v] = ...
  capture_hyp(goal_id, orbit, arr_time, park_r, origin_coe, v_in)
% CAPTURE_HYP(goal_id, orbit, arr_time, park_r, park_i, origin_coe, v_in)
%   computes the trajectory the spacecraft will follow
%   to enter the sphere of influence of object GOAL_ID at time
%   ARR_TIME and to position itself into a circular parking orbit of 
%   radius PARK_R.
%
%   It uses the last points of the arrival interplanetary orbit 
%   (computed via the patched conics method) stored in ORBIT to generate
%   a trajectory with orbital elements ORIGIN_COE that will allow it
%   to reach velocity V_IN at the entrance of the body's SOI.
%
%   [traj, delta_v] = CAPTURE_HYP(...) returns the computed capture
%   trajectory TRAJ and the needed change of velocity DELTA_V.
%
%   Options for the Sun are not contemplated since that would be
%   the general case of an interplanetary trajectory.
%
%   goal_id  - identifier of the destination planet:
%                    1 = Mercury
%                    2 = Venus
%                    3 = Earth
%                    4 = Mars
%                    5 = Jupiter
%                    6 = Saturn
%                    7 = Uranus
%                    8 = Neptune
%                    9 = Pluto
%                   10 = Vesta
%                   11 = Ceres
%
%   orbit    - last two points of the interplanetary trajectory
%              computed via the patched conics method
%
%   arr_time - array specifying time of arrival with elements 
%                    (in this order):
%                     year         - range: 1901 - 2099
%                     month        - range: 1 - 12
%                     day          - range: 1 - 31
%                     hour         - range: 0 - 23
%                     minute       - range: 0 - 60
%                     second       - range: 0 - 60
%
%   park_r   - radius of the circular parking orbit around the body
%
%   origin_coe - classical orbital elements of the origin interplanetary
%                 orbit:
%                h    = angular momentum (km^2/s)
%                e    = eccentricity
%                RA   = right ascension of the ascending
%                       node (rad)
%                incl = inclination of the orbit (rad)
%                w    = argument of perigee (rad)
%                TA   = true anomaly (rad)
%                a    = semimajor axis (km)
%
%   v_in - velocity of the spacecraft at the entry into the SOI
%          of the planet

    %% Argument validation
    validateattributes(arr_time,{'double'},{'size',[1 6]})
    validateattributes(origin_coe,{'double'},{'size',[1 7]})

    %% Data
    global mu

    masses = 10^24 * [0.330104
                      4.86732
                      5.97219
                      0.641693
                      1898.13
                      568.319
                      86.8103
                      102.410
                      0.01309
                      0.000259
                      0.0009393
                      1989100]; %[kg]

	radii = [2439.7
             6051.8 
             6371
             3389.5
             69911
             58232
             25362
             24622
             1151
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
                 353649000000
                 413690250];%[km]
    
    G = 6.6742e-20; %[km^3/kg/s^2]
    
    %% Input data
    %SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
    pl_SOI = (masses(goal_id)/masses(12))^(2/5)...
        * distances(goal_id); %[km]
    
    pl_mu = G * masses(goal_id); %[km^3/s^2]
    pl_radius = radii(goal_id); %[km]
    
    %% Needed variables
    [~, pl_r0, v_arr, ~] =...
        planet_elements_and_sv(goal_id, arr_time(1),arr_time(2),...
                        arr_time(3),arr_time(4),arr_time(5),arr_time(6));

    vinf = norm(v_arr-v_in); % [km/s^2] , v-infinity

    rp = pl_radius + park_r; % [km] perigee
    e = 1 + rp*vinf^2/pl_mu; % eccentricity
    a = rp/(e-1); % [km] semi-major axis
    
    % target radius for the right hyperbola (Fig. 8.13-8.14 Curtis)
    Delta = rp*sqrt(1+2*pl_mu/(rp*vinf^2)); % [km] aiming radius
    
    %Angle between arrival and departure branch of the hyperbola
    half_delta = asin(1/e);
    
    v_hyp = sqrt(vinf^2+2*pl_mu/rp);%sqrt(-2*pl_mu/rp + pl_mu/a);
    %vc = sqrt(2*pl_mu/rp);? Eq 8.59 with circular parking orbit
    vc = sqrt(pl_mu/rp); % velocity of parking orbit

    h = Delta*vinf; % [km^2/s] specific angular momentum
    RA = origin_coe(3); %[rad] right ascension of ascending node
    incl = origin_coe(4); %[rad] inclination

    n = sqrt(pl_mu/a^3); % mean motion
    
    %% Trajectory computation
    hyp = zeros(100*24*3600/60,3);
    
    %Angle and direction of the orbit at the entering point
    in_dir = (orbit(1,1:3)-orbit(2,1:3))'; %exit vector: (1,1:3)<-(2,1:3)
    in_angle = deg2rad(atan2d_0_360(in_dir(2),in_dir(1)));

    %Initial point
    coe = [h, e, RA, incl, 0, 0];
    [rprova,~] = sv_from_coe(coe, pl_mu);
    alpha = deg2rad(atan2d_0_360(rprova(2),rprova(1)));
    
    %Desired characteristics to align with the interplanetary orbit
    xi_des = in_angle - pi - half_delta;
    alpha_des = pi/2 + xi_des;
    w_des = alpha_des - alpha;        
    
    counter = 1;
    for t=0:60:100*24*3600
        %Using Kepler's method to compute the point
        M = n*t;
        F = kepler_H(e,M);
        cosf = (cosh(F)-e)/(1-e*cosh(F)); %(Eq. 3.41b) Curtis
        f = acos(cosf);
        
        if (goal_id == 11) %to ensure desired alignment
            coe = [h, e, RA, incl, w_des+pi,-f];
        else 
            coe = [h, e, RA, incl, w_des+pi,f];
        end
        [r,~] = sv_from_coe(coe, pl_mu);
        
        %Sometimes the above algorithm produces NaN elements because of the
        %kepler_H function (it seems not to be able to deal with high
        %numbers)
        if(any(isnan(r)))
            if(hyp(2,:) ~= [0 0 0]) %from the third point onward
                diff = hyp(counter-1,:)-hyp(counter-2,:);
                diff = 60*norm(v_in)*diff/norm(diff);
                point = hyp(counter-1,:)' + diff';
            else %first two points
                if (goal_id > 4) %for outer planets/elements
                    coe = [h, e, RA, incl, w_des+pi,-t/6];
                else %for inner planets/elements
                    coe = [h, e, RA, incl, w_des,t/6];
                end
                [peri,~] = sv_from_coe(coe, pl_mu);
                point = pl_r0' + peri';
            end
        else %if kepler_H returned a valid result
             point = pl_r0' + r';
        end
        
        %Adding the point to the list
        hyp(counter,:) = point';
        counter = counter+1;
        
        %To stop computing when the SOI is exited
        if (all(hyp(1,:) ~= [0 0 0]) && norm(point'-hyp(1,:))>= pl_SOI)
            break;
        end
    end
    %Getting rid of unused elements in the array
    hyp = hyp(1:counter-2, 1:3);
    
    %% Hyperbola plot
    plot3(hyp(:,1),hyp(:,2),hyp(:,3),'m-')
    
    %% Output arguments
    traj = flip(hyp);
    delta_v = v_hyp - vc;
end    