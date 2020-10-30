function orb = intpl_orbit2(arr_days, r_init, v_init)
%   INTERPLANETARY_ORBIT(arr_days, r_init, v_init) computes the 
%   interplanetary orbit a spacecraft of 1000kg would travel
%   around the Sun in ARR_DAYS days, starting at position R_INIT
%   with velocity V_INIT.
%   It uses the patched conics method.
%
%   orb = INTPL_ORBIT(...) returns the points composing the orbit.
%
%   It uses rkf45 to numerically integrate Equation 2.22 in
%   "Orbital Mechanics for Engineering Students" - Howard D. Curtis.
%
%   arr_days  - days needed for the spacecraft to arrive at destination
%
%   r_init    - initial position of the spacecraft
%
%   v_init    - initial velocity of the spacecraft
%
% User M-function required:   rkf45
% User subfunctions required: rates, output

    %% Constants
    global mu
    hours = 3600; %[s]  

    %% Input data:
    
    r0 = r_init; %[km,km,km]
    v0 = v_init; %[km/s,km/s,km/s]

    t0 = 0; %[s]
    tf = arr_days*24*hours; % [s]

    Sun_mu  = mu; %[km^3/s^2]
    
    %% Numerical integration:
    
    %Initial condition: position, velocity
    y0    = [r0 v0]';
    [t,y] = rkf45(@rates, [t0 tf], y0, 1.e-20);

    %% Output the results:
    output
    
    orb = y;

    return 

    %% ~~~~~~~~~~~~~~~~~~~~~~~~
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
        
        x    = f(1);
        y    = f(2);
        z    = f(3);
        vx   = f(4);
        vy   = f(5);
        vz   = f(6);

        r    = norm([x y z]);

        ax   = -Sun_mu*x/r^3;
        ay   = -Sun_mu*y/r^3;
        az   = -Sun_mu*z/r^3;

        dydt = [vx vy vz ax ay az]';    
    end %rates
    % ~~~~~~~~~~~~~~~~~~~~~~~


    %% ~~~~~~~~~~~~~~~~~~~~~~~
    function output
    %{
      This function computes the maximum and minimum radii, the times they
      occur and and the speed at those times. It prints those results to
      the command window and plots the orbit.
    %}
        
%         %% Output to the command window:
%         fprintf('\n\n--------------------------------------------------------\n')
%         fprintf('\n Spacecraft Orbit\n')
%         fprintf('\n The initial position is [%g, %g, %g] (km).',...
%                                                              r0(1), r0(2), r0(3))
%         fprintf('\n   Magnitude = %g km\n', norm(r0))
%         fprintf('\n The initial velocity is [%g, %g, %g] (km/s).',...
%                                                              v0(1), v0(2), v0(3))
%         fprintf('\n   Magnitude = %g km/s\n', norm(v0))
%         fprintf('\n--------------------------------------------------------\n\n')
% 
%         %% Plot the results:
%         
%         hold on
%         plot3(  y(:,1),    y(:,2),    y(:,3), 'k', 'LineWidth', 2)

    end %output
    % ~~~~~~~~~~~~~~~~~~~~~~~

end %orbit