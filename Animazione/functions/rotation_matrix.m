function [R_x,R_y,R_z] = rotation_matrix(alpha_x,beta_y,gamma_z)
    deg = pi/180;
    R_x = [1 0 0; 0 cos(alpha_x*deg) sin(-alpha_x*deg); 0 sin(alpha_x*deg) cos(alpha_x*deg);];
    R_y = [cos(beta_y*deg) 0 sin(beta_y*deg); 0 1 0; -sin(beta_y*deg) 0 cos(beta_y*deg)];
    R_z = [cos(gamma_z*deg) -sin(gamma_z*deg) 0; sin(gamma_z*deg) cos(gamma_z*deg) 0; 0 0 1];
end

