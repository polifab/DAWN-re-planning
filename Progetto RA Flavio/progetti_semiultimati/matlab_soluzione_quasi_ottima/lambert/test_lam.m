r1 = [terra.pos(11483,1); terra.pos(11483,2); 0];
r2 = [venere.pos(21843,1);venere.pos(21843,2);0];

u = (6.67e-11*1.9891e30)*1e-9;

[V1, V2, extremal_distances, exitflag] = lambert(r1', r2', tf, 1, u);