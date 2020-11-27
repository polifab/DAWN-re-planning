function [orb, t, deltav] = ...
                    park_orbit_change(obj_id, body_pos, r1, r2, tf, grade)
% PARK_ORBIT_CHANGE computes the change of parking orbit around a body.
%
%   INPUT:
%       obj_id   - identifier of the considered object
%       body_pos - position of the object with respect to the Sun
%       r1       - spacecraft position in the departure parking orbit
%       r2       - spacecraft position in the arrival parking orbit
%       tf       - time of flight
%       grade    - string to set retrograde or prograde orbit
%   OUT:
%       orb      - orbit points
%       t        - time istants for each point

    %% Data
    
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
    
    G = 6.6742e-20; %[km^3/kg/s^2]
    
    %% Setting problem

    r2 = r2 - body_pos; % arrival position
    r1 = r1 - body_pos; % departure position

    mu = G*masses(obj_id);

    [V1, V2] = lambert(r1, r2, tf, grade, mu);

    %% Integration
    
    y0 = [r1 V1]';
    [t,y] = rkf45(@rates, [0 tf], y0);

    %% Plotting
    
    new_r1 = body_pos + y(1,1:3);
    new_r2 = body_pos + y(end,1:3);

    plot3(new_r1(:,1), new_r1(:,2), new_r1(:,3), 'o')
    plot3(new_r2(:,1), new_r2(:,2), new_r2(:,3), '*')

    %% Delta-v computation

    R2 = sqrt(dot(r2,r2)) ;           % module of r2
    R1 = sqrt(dot(r1,r1)) ;           % module of r1
    v_initial = sqrt(dot(V1,V1));     % module of V1
    v_final = sqrt(dot(V2,V2));       % module of V2

    v_park1 = sqrt(mu/R1);             % initial circular orbit velocity 
    v_park2 = sqrt(mu/R2);             % final circular orbit velocity

    deltav1 = v_initial - v_park1;     % initial delta v
    deltav2 = v_final - v_park2;       % final delta v


    %% Setting output parameters
    
    orb  = body_pos +  y(:,1:3);
    plot3(orb(:,1), orb(:,2), orb(:,3), '-.', 'LineWidth', 2)
    deltav = abs(deltav1) + abs(deltav2);   % total delta v
    
    %% ~~~~~~~~~~~~~~~~~~~~~~~
    function dydt = rates(t,f)
    %{
      This function calculates the acceleration vector using Equation 2.22

      t          - time
      f          - column vector containing the position vector and the
                   velocity vector at time t
      x, y, z    - components of the position vector r
      r          - the magnitude of the the position vector
      vx, vy, vz - components of the velocity vector v
      ax, ay, az - components of the acceleration vector a
      dydt       - column vector containing the velocity and acceleration
                   components
    %}
        mu_vesta = 17.8;
        x    = f(1);
        y    = f(2);
        z    = f(3);
        vx   = f(4);
        vy   = f(5);
        vz   = f(6);

        r    = norm([x y z]);

        ax   = -mu_vesta*x/r^3;
        ay   = -mu_vesta*y/r^3;
        az   = -mu_vesta*z/r^3;

        dydt = [vx vy vz ax ay az]';
        
    end %rates
    % ~~~~~~~~~~~~~~~~~~~~~~~
end