Vesta_r3 = 1.0e+08 * [1.964223506918297 -2.676863670878835 -0.158843745137528];

r2 = (1.0e+08 * [ 1.964223326309674  -2.676868093233760  -0.158843836939114]) - Vesta_r3;
r1 = (1.0e+08 * [ 1.964223456944477  -2.676854348001822  -0.158844017768201]) - Vesta_r3;


%r1 = 1.0e+02 * [8.459989460540838   3.757532311395191  -1.141148941063138];
%r2 = 1.0e+02 * [-1.004253938250492   4.298513159856268   0.335401480537469];

tf = 15000;
global mu
mu = 17.8;

[V1, V2] = lambert(r1, r2, tf, 'pro');
%[V1, V2, extremal_distances, exitflag] = lambert_rodyo(r1, r2, tf, 0, mu);
%V1_new = [V1(2) V1(1) V1(3)];
y0 = [r1 V1]';
[t,y] = rkf45(@rates, [0 tf], y0);
%%
new_r1 = Vesta_r3 + r1;
new_r2 = Vesta_r3 + r2;

orb  = Vesta_r3 +  y(:,1:3);

plot3(new_r1(:,1), new_r1(:,2), new_r1(:,3), 'o')
hold on;
plot3(new_r2(:,1), new_r2(:,2), new_r2(:,3), '*')
plot3(orb(:,1), orb(:,2), orb(:,3), '-.', 'LineWidth', 2)


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
