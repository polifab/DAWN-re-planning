function [traj, delta_v] = escape_hyp(obj_id, orbit,...
                               dep_time, park_r, goal_coe, v_out)
% ESCAPE_HYP(planet_id,goal_id,orbit,dep_time,park_r,park_i,goal_coe,v_out)
%   computes the trajectory the spacecraft will follow
%   to escape the sphere of influece of the object OBJ_ID AT
%   DEP_TIME.
%
%   It uses the first points of the departing interplanetary orbit 
%   (computed via the patched conics method) stored in ORBIT and the
%   parking orbit radius PARK_R and inclination PARK_I to generate
%   a trajectory with orbital elements GOAL_COE that will allow it
%   to reach velocity V_OUT at the end of the travel.
%
%   [traj, delta_v] = ESCAPE_HYP(...) returns the computed escape
%   trajectory TRAJ and the needed change of velocity DELTA_V.
%
%   Options for the Sun are not contemplated since that would be
%   the general case of an interplanetary trajectory.
%
%   obj_id  - identifier of the origin planet:
%                            1 = Mercury
%                            2 = Venus
%                            3 = Earth
%                            4 = Mars
%                            5 = Jupiter
%                            6 = Saturn
%                            7 = Uranus
%                            8 = Neptune
%                            9 = Pluto
%                           10 = Vesta
%                           11 = Ceres
%
%   orbit    - first two points of the interplanetary trajectory
%              computed via the patched conics method
%
%   dep_time - array specifying time of departure with elements 
%                    (in this order):
%                     year         - range: 1901 - 2099
%                     month        - range: 1 - 12
%                     day          - range: 1 - 31
%                     hour         - range: 0 - 23
%                     minute       - range: 0 - 60
%                     second       - range: 0 - 60
%
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
%
%   v_out - escape velocity from the object SOI

    %% Argument validation
    validateattributes(dep_time,{'double'},{'size',[1 6]})
    validateattributes(goal_coe,{'double'},{'size',[1 7]})

    %% Constants
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
    pl_SOI = (masses(obj_id)/masses(12))^(2/5)...
        * distances(obj_id); %[km]
    
    pl_mu = G * masses(obj_id); %[km^3/s^2]
    pl_radius = radii(obj_id); %[km]
    
    park_i = goal_coe(4);

    %% Needed variables
    [~, pl_r0, v_dep, ~] =...
        planet_elements_and_sv(obj_id,dep_time(1),dep_time(2),...
                        dep_time(3),dep_time(4),dep_time(5),dep_time(6));
   
%     V_dep = norm(v_dep); %[km/s], norm of departing velocity

    %v-infinity of the departure hyperbola
    vinf = norm(v_out - v_dep); %[km/s]

    %Hyperbola characteristics
    rp = pl_radius + park_r; %[km], periapsis
    e = 1 + rp*vinf^2/pl_mu; %eccentricity
    a = rp/(e-1); %[km], semi-major axis
%     b = a*sqrt(e^2-1); %[km], semi-minor axis

    %Velocity at hyperbola periapsis
    vp = sqrt(vinf^2 + 2*pl_mu/rp);

    %Angle between hyperbola center and exiting-branch
%     beta = acos(1/e);
    half_delta = asin(1/e);

    %Hyperbola orbital elements
    h = rp*vp;
    RA = goal_coe(3);
    incl = goal_coe(4);
%     w = goal_coe(5);

    %Mean motion
    n = sqrt(pl_mu/a^3);
    
    %Gravitational parameter of the origin planet
    mu_dep = masses(obj_id) * G; %[km^3/s^2]
    %Velocity of the spacecraft after burn
    v_b = sqrt(vinf^2 + 2*mu_dep/park_r); %[km/s]
    %Velocity of circular parking orbit around origin planet
    v_park = sqrt(mu_dep/park_r);
    
    %% Trajectory computation
    rr = [];
    out_dir = (orbit(2,1:3)-orbit(1,1:3))'; %exit vector: (2,1:3)<-(1,1:3)
    out_angle = deg2rad(atan2d_0_360(out_dir(2),out_dir(1)));
%     up_angle = deg2rad(-atan2d_0_360(out_dir(3),-out_dir(1)));

%     p = -a*(1-e^2); %semilatum [km]
%     tran = acos(1/e-rp/p); %true anomaly at perigee
    
    coe = [h, e, RA, incl, 0, 0];
    [rprova,~] = sv_from_coe(coe, pl_mu);
    alpha = deg2rad(atan2d_0_360(rprova(2),rprova(1)));

    xi_des = out_angle - pi - half_delta;
    alpha_des = pi/2 + xi_des;
    w_des = alpha_des - alpha;        
    
    for t=0:60:100*24*3600%24*3600
        M = n*t; %Hyperbolic mean anomaly
        F = kepler_H(e,M); %Hyperbolic eccentric anomaly
        cosf = (e-cosh(F))/(e*cosh(F)-1);
        f = acos(cosf); %True anomaly
        coe = [h, e, RA, incl, w_des, f];
        [r,~] = sv_from_coe(coe, pl_mu); %spacecraft position
        if(size(rr,1)>1 && any(isnan(r)))
            diff = rr(end,:)-rr(end-1,:);
            point = rr(end,:)' + diff';
        else
             point = pl_r0' + r';
        end
        rr = cat(1,rr,point');
        if size(rr,1)>0 && norm(point'-rr(1,:))>= pl_SOI
            break;
        end
    end

    %% Deleted because useless, I think
    %Angle of orientation of escape trajectory, to be aligned with
    %the escape velocity vector
%     out_dir = Rotz(goal_coe(3))'*Rotx(park_i)'*...
%         (orbit(2,1:3)-orbit(1,1:3))'; %exit vector: (2,1:3)<-(1,1:3)
%     out_angle = deg2rad(atan2d_0_360(out_dir(2),out_dir(1)));
% 
%     t = 0:0.1:5;
% 
% %     Parametric hyperbola equations
%     xh_l = -a*cosh(t);
%     xh_r = a*cosh(t);
%     yh = b*sinh(t);
% 
%     hyp = [];
%     for i = 1:length(t)
%         point = pl_r0' + Rotz(goal_coe(3))*Rotx(park_i)*...
%                     Rotz(out_angle)*Rotz(beta)*...
%                     ([xh_l(i); -yh(i);0] + [-(rp-a);0;0]);
%         hyp = cat(1,hyp,point');
%         if norm(hyp(size(hyp,1),:)-hyp(1,:))>= pl_SOI
%             break;
%         end
%     end
    
    %% Hyperbola plot
    plot3(rr(:,1),rr(:,2),rr(:,3),'m-')%plot3(hyp(:,1),hyp(:,2),hyp(:,3),'m-')
    
    %% Output arguments
    traj = rr;%hyp;
    delta_v = v_b - v_park;
end