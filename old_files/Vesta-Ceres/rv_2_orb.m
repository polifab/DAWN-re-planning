% function that computes the classical orbital elements [a e i RAAN w f] 
% from the state vector (r,v) referred to a body with standard gravitational 
% parameters mu

function [a e i RAAN w f] = rv_2_orb(r,v,mu)

r0 = sqrt(dot(r,r));    %module of r
v0 = sqrt(dot(v,v));    %module of v

E = 0.5*v0^2 - mu/r0;   %specific orbital energy
a = -mu/(2*E);          % 1 - semi major axis

h = cross(r,v);         %specific angular momentum vector
h0 = sqrt(dot(h,h));    %module of h
e_v = (cross(v,h)/mu-r/r0); %eccentricity vector

e = sqrt(dot(e_v,e_v));       % 2 - eccentricity  (e = sqrt(1+(2*E*h0^2/mu^2)))

i = acos(dot([0 0 1]',h)/h0)  % 3 - inclination

n_vett = cross([0 0 1]', h);      %line of nodes vector
n_mod = sqrt(dot(n_vett,n_vett)); % module of n_vett
n = n_vett/n_mod;                 %line of nodes versor

% 4 - RAAN 
if dot(n,[0 1 0]) > 0
    RAAN = acos(dot([ 1 0 0],n))
    
else if dot(n,[0 1 0]) < 0             % case dot(n,[0 1 0])= 0 is not considered
    RAAN = 2*pi - acos(dot([1 0 0],n)) 
    end 
end

% 5 - w
 if (dot(e_v,[0 0 1]) > 0)
     w = acos(dot(n,e_v)/e)
     
 else if (dot(e_v,[0 0 1]) < 0)             % case dot(e_v,[0 0 1])= 0 is not considered
     w = 2*pi - acos(dot(n,e_v)/e)
     end
 end
 
% 6 - true anomaly
if (dot(r,v) > 0)
    f = acos(dot(e_v,r)/(e*r0))
else if (dot(r,v) < 0)                    % case dot(r,v)= 0 is not considered
    f = 2*pi - acos(dot(e_v,r)/(e*r0))
    end

end