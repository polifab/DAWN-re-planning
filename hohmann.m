function [deltav, deltat] = hohmann(r1, r2, mu)
% HOHMANN return the time window and velocity change needed to change
%   parking orbit around a body, following a Hohmann transfer.
%
%   r1 - radius of the first parking orbit
%   r2 - radius of the second parking orbit
%   mu - gravitational parameter of the body
%

    vinit = sqrt(mu/r1);                         % first circular orbit velocity
    vtrans_p = sqrt(2*mu*r2/(r1*(r1 + r2)));     % initial velocity transfer (perielion)
    deltav1 = vtrans_p - vinit;                  % initial delta v
    
    vfin = sqrt(mu/r2);                          % final circular orbit velocity
    vtrans_a = sqrt(2*mu*r1/(r2*(r1 + r2)));     % final velocity transfer (aphelion)
    deltav2 = vtrans_a - vfin;                   % final delta v
    
    deltav = deltav1 + deltav2;                  % total delta v
    deltat = pi*sqrt((r1 + r2)^(3)/(8*mu));      % total delta t
end