function [eanom, tanom] = kepler2 (manom, ecc)

% solve Kepler's equation for circular,
% elliptic and hyperbolic orbits

% Danby's method with Mikkola's initial guess

% input

%  manom = mean anomaly (radians)
%  ecc   = orbital eccentricity (non-dimensional)

% output

%  eanom = eccentric anomaly (radians)
%  tanom = true anomaly (radians)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global niter

% define convergence criterion

ktol = 1.0e-10;

pi2 = 2 * pi;

if (ecc == 0)
    
    % circular orbit

    tanom = manom;
    
    eanom = manom;
    
    return;
    
end

den = 1 / (4 * ecc + 0.5);

if (ecc < 1)
    
    xma = manom - pi2 * fix(manom / pi2);

    alpha = (1 - ecc) * den;
    
else
    
    xma = manom;
    
    alpha = (ecc - 1) * den;
    
end

beta = 0.5 * xma * den;

z = (sqrt(alpha * alpha * alpha + beta * beta) + beta) ^ (1/3);

s = 2 * beta /(z * z + alpha + alpha * alpha / (z * z));

if (ecc > 1)
    
    % hyperbolic orbit
    
    ds = 0.071 * s ^ 5 / (ecc * (1 + 0.45 * s ^ 2) * (1 + 4 * s ^ 2));
    
else
    
    % elliptic orbit
    
    ds = -0.078 * s * s * s * s * s / (1 + ecc);
    
end

s = s + ds;

% initial guess

if (ecc > 1)
    
    % hyperbolic orbit
    
    eanom = 3 * log(s + sqrt(1 + s ^ 2));
    
else
    
    % elliptic orbit
    
    eanom = xma + ecc * (3 * s - 4 * s * s * s);
    
end

% perform iterations

niter = 0;

while(1)
    
    if (ecc < 1)

        % elliptic orbit

        s = ecc * sin(eanom);
        c = ecc * cos(eanom);

        f = eanom - s - xma;
        fp = 1 - c;
        fpp = s;
        fppp = c;
        
    else

        % hyperbolic orbit

        s = ecc * sinh(eanom);
        c = ecc * cosh(eanom);

        f = s - eanom - xma;
        fp = c - 1;
        fpp = s;
        fppp = c;
        
    end

    niter = niter + 1;

    % check for convergence

    if (abs(f) <= ktol || niter > 20)
        
        break;
        
    end

    % update eccentric anomaly

    delta = -f / fp;

    deltastar = -f / (fp + 0.5 * delta * fpp);

    deltak = -f / (fp + 0.5 * deltastar * fpp ...
        + deltastar * deltastar * fppp / 6);

    eanom = eanom + deltak;
    
end

if (niter > 20)
    
    clc; home;
    
    fprintf('\n\n   more than 20 iterations in kepler2 \n\n');
    
    keycheck;
    
end

% compute true anomaly

if (ecc < 1)
    
    % elliptic orbit

    sta = sqrt(1 - ecc * ecc) * sin(eanom);
    
    cta = cos(eanom) - ecc;
    
else
    
    % hyperbolic orbit

    sta = sqrt(ecc * ecc - 1) * sinh(eanom);
    
    cta = ecc - cosh(eanom);
    
end

tanom = atan3(sta, cta);


