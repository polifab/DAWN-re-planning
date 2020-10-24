function [traj, delta_v] = ...
  capture_hyp(goal_id, orbit, arr_time, park_r, origin_coe, v_in)
% ESCAPE_HYP(goal_id, orbit, arr_time, park_r, park_i, origin_coe, v_in)
%   computes the trajectory the spacecraft will follow
%   to enter the sphere of influence of object GOAL_ID at time
%   ARR_TIME and to position itself into a circular parking orbit
%   of radius PARK_R and inclination PARK_I.
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
    
    park_i = origin_coe(4);

    %% Needed variables
    [~, pl_r0, v_arr, ~] =...
        planet_elements_and_sv(goal_id, arr_time(1),arr_time(2),...
                        arr_time(3),arr_time(4),arr_time(5),arr_time(6));

    vinf = norm(v_arr-v_in);%norm(v_in-v_arr);

    rp = pl_radius + park_r;
    e = 1 + rp*vinf^2/pl_mu;
    a = rp/(e-1);
    b = a*sqrt(e^2-1);
    
    % target radius for the right hyperbola
    Delta = rp*sqrt(1+2*pl_mu/(rp*vinf^2)); %Delta = sqrt(a*(1 - e^2)*pl_mu/vinf);
    
    %Angle between arrival and departure branch of the hyperbola
    half_delta = asin(1/e);
    
    v_hyp = sqrt(vinf^2+2*pl_mu/rp);%sqrt(-2*pl_mu/rp + pl_mu/a);
    %vc = sqrt(2*pl_mu/rp);? Eq 8.59 with circular parking orbit
    vc = sqrt(pl_mu/rp); % velocity of parking orbit

    beta = acos(1/e);
    
    h = Delta*vinf;
    RA = origin_coe(3); %[rad]
    incl = origin_coe(4); %[rad]
    w = origin_coe(5); %[rad]
    
    TA = origin_coe(6); %[rad]

    n = sqrt(pl_mu/a^3);
    
    %% Trajectory computation
    rr = zeros(100*24*3600/60,3);
    
    in_dir = (orbit(1,1:3)-orbit(2,1:3))'; %exit vector: (1,1:3)<-(2,1:3)
    in_angle = deg2rad(atan2d_0_360(in_dir(2),in_dir(1)));

    coe = [h, e, RA, incl, 0, 0];
    [rprova,~] = sv_from_coe(coe, pl_mu);
    alpha = deg2rad(atan2d_0_360(rprova(2),rprova(1)));
    
    xi_des = in_angle - pi - half_delta;
    alpha_des = pi/2 + xi_des;
    w_des = alpha_des - alpha;        
    counter = 1;
    for t=0:60:100*24*3600%ceil(pl_SOI/norm(v_in))
        M = n*t;
        F = kepler_H(e,M);
        cosf = (cosh(F)-e)/(1-e*cosh(F));%Eq 3.41b %cosf = (e-cosh(F))/(e*cosh(F)-1);
        f = acos(cosf);
        coe = [h, e, RA, incl, w_des, f];
        [r,~] = sv_from_coe(coe, pl_mu);
        if(any(isnan(r)))
            if(rr(2,:) ~= [0 0 0])
                diff = rr(end,:)-rr(end-1,:);
                point = rr(end,:)' + diff';
            else
%                 diff = Rotz(RA)*Rotx(incl)*[round(norm(v_in)*t);0;0];
%                 point = pl_r0' + rprova' + diff';
%                 t = t + 24*3600;
                coe = [h, e, RA, incl, w_des, t/6];
                [peri,~] = sv_from_coe(coe, pl_mu);
                point = pl_r0' + peri';
            end
        else
             point = pl_r0' + r';
        end
        rr(counter,:) = point';
        counter = counter+1;
        if all(rr(1,:) ~= [0 0 0]) && norm(point'-rr(1,:))>= pl_SOI
            break;
        end
    end
    rr = rr(1:counter-2, 1:3);

    %% Deleted because useless, I think
    %Angle of orientation of escape trajectory
%     in_dir = Rotz(origin_coe(3))'*Rotx(park_i)'*...
%         (orbit(1,1:3)-orbit(2,1:3))'; %entry vector: (end-1,1:3)<-(end,1:3)
%     in_angle = deg2rad(atan2d_0_360(in_dir(2),in_dir(1)));
% 
%     t = 0:0.1:5;
% 
%     %Parametric hyperbola equations
%     xh_l = -a*cosh(t);
%     xh_r = a*cosh(t);
%     yh = b*sinh(t);
%     
%     hyp = [];
%     for i = 1:length(t)
%         point = pl_r0' + Rotz(origin_coe(3))*Rotx(park_i)*...
%                     Rotz(in_angle)*Rotz(2*half_delta+beta)*...
%                     ([xh_l(i); -yh(i);0] + [-(rp-a);0;0]);
%         hyp = cat(1,hyp,point');
%         if norm(hyp(size(hyp,1),:)-hyp(1,:))>= pl_SOI
%             break;
%         end
%     end
    
    %% Hyperbola plot
    plot3(rr(:,1),rr(:,2),rr(:,3),'b-')%plot3(hyp(:,1),hyp(:,2),hyp(:,3),'b-')
    
    %% Output arguments
    traj = rr;%flip(hyp);
    delta_v = v_hyp - vc;
end    