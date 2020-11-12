function [orb, t, deltav] = cambio_orbita_park(Vesta_r3, r1, r2, tf, grade)
%% cambio_orbita_park
%  IN:
%  Cambio di orbita di parcheggio su Vesta (da adattare per generalizzare)
%  Vesta_r3: posizione di vesta rispetto al Sole in quel momento
%  r1: posizione di partenza dell'orbita di parcheggio
%  r2: posizione di arrrivo
%  tf: tempo di volo
%  grade: parametro per un'orbita prograda o retrograda per il cambio di orbita
%  OUT:
%  orb: orbita del cambio di orbita di parcheggio
%  t: tempi in cui è valutata l'orbita
%% Setting problem

    %Vesta_r3 = 1.0e+08 * [1.964223506918297 -2.676863670878835 -0.158843745137528]; % Vesta position respect to Sun at arrival

    r2 = r2 - Vesta_r3; % arrival pos
    r1 = r1 - Vesta_r3; % dep pos

    %tf = 12000;
    global mu
    mu = 17.8;

    [V1, V2] = lambert(r1, r2, tf, grade);

    y0 = [r1 V1]';
    [t,y] = rkf45(@rates, [0 tf], y0);

 
%% plotting
    new_r1 = Vesta_r3 + r1;
    new_r2 = Vesta_r3 + r2;

    orb  = Vesta_r3 +  y(:,1:3);

    plot3(new_r1(:,1), new_r1(:,2), new_r1(:,3), 'o')
    hold on;
    plot3(new_r2(:,1), new_r2(:,2), new_r2(:,3), '*')
    plot3(orb(:,1), orb(:,2), orb(:,3), '-.', 'LineWidth', 2)

 %% deltav computation
 
R2 = sqrt(dot(r2,r2)) ;           % module of r2
R1 = sqrt(dot(r1,r1)) ;           % module of r1
v_initial = sqrt(dot(V1,V1));     % module of V1
v_final = sqrt(dot(V2,V2));       % module of V2

v_park1 = sqrt(mu/R1)             % initial circular orbit velocity 
v_park2 = sqrt(mu/R2)             % final circular orbit velocity

deltav1 = v_initial - v_park1     % initial delta v
deltav2 = v_final - v_park2       % final delta v

deltav = abs(deltav1) + abs(deltav2)   % total delta v

    
    mu = 1.327565122000000e+11; %[km^3/s^2]
    
end
%%
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
        mu_vesta = 17.8; %1.327565122000000e+11
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
