function [inclination_change, true_anomaly, Dv] = change_of_plane(i1, i2, W1, W2, v_initial)
%   Function for computing relevant quantities in change of plane
%
%   INPUT:
%       i1 : initial orbital inclination (rad)  
%       i2 : target orbital inclination  (rad)
%       W1 : initial Right Ascension of the Ascending Node (RAAN) (rad)
%       W2 : target RAAN (rad)
%       v_initial : velocity in initial orbit (km/s)s
%
%   OUTPUT:
%       inclination_change  : change of inclination between initial and target
%       true_anomaly        : position of spacecraft in orbit
%       DV                  : delta v required for change of plane
%
%


    theta   = acos(cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(W2 - W1));

    inclination_change  = theta;
    
    true_anomaly        = acos((cos(i1)*cos(theta) - cos(i2))/(sin(i1)*sin(theta)));
    
    Dv      = 2*v_initial*sin(theta/2);



end
