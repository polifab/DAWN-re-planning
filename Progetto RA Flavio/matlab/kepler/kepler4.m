function [r, tanom] = kepler4(t, q, ecc)

% Kepler's equation for heliocentric
% parabolic and near-parabolic orbits

% input

%  t   = time relative to perihelion passage (days)
%  q   = perihelion radius (AU)
%  ecc = orbital eccentricity (non-dimensional)

% output

%  r     = heliocentric distance (AU)
%  tanom = true anomaly (radians)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kgauss = 0.01720209895;

e2 = 0;

fac = 0.5 * ecc;

i = 0;

tau = kgauss * t;

while(1)
    
    i = i + 1;

    e20 = e2;

    a = 1.5 * sqrt(fac / (q * q * q)) * tau;

    b = (sqrt(a * a + 1) + a)^(1/3);

    u = b - 1 / b;

    u2 = u * u;

    e2 = u2 * (1 - ecc) / fac;

    [c1, c2, c3] = stumpff(e2);

    fac = 3 * ecc * c3;

    % check for convergence

    if (abs(e2 - e20) < 1e-9) || (i > 15)
        
        break;
        
    end
    
end

if (i > 15)
    
    clc; home;
    
    fprintf('\n\n   more than 15 iterations in kepler4 \n\n');
    
    keycheck;
    
end

% heliocentric distance (AU)

r = q * (1 + u2 * c2 * ecc / fac);

% x and y coordinates

x = q * (1 - u2 * c2 / fac);

y = q * sqrt((1 + ecc) / fac) * u * c1;

% true anomaly (radians)

tanom = atan3(y, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c1, c2, c3] = stumpff(e2)

% Stumpff functions for argument E^2

% input

%  e2 = eccentric anomaly squared

% output

%  c1, c2, c3 = Stumpff functions

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1 = 0;
c2 = 0;
c3 = 0;

deltac = 1;

n = 1;

while (1)
    
    c1 = c1 + deltac;

    deltac = deltac / (2 * n);

    c2 = c2 + deltac;

    deltac = deltac / (2 * n + 1);

    c3 = c3 + deltac;

    deltac = -e2 * deltac;

    n = n + 1;

    % check for convergence

    if (abs(deltac) < 1e-12)
        
        break;
        
    end
    
end


