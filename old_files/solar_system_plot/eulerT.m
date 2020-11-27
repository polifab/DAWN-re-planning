function [M_rot] = eulerT(W,w,i)

M1 = [cos(W)    -sin(W) 0;
      sin(W)    cos(W)  0;
      0         0       1];
M2 = [1         0       0;
      0         cos(i)  -sin(i);
      0         sin(i)  cos(i)];
M3 = [cos(w)    -sin(w) 0;
      sin(w)    cos(w)  0;
      0         0       1];
  
M_rot = M1*M2*M3;
end