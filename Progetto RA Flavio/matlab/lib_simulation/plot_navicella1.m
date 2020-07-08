%conosco la struttura vol1 che contiene la posizione
vol1.pnav = eye(3,1);
vol1.manom = mean_anom_t(vol1.mean_anom0, venere.u_sol ,vol1.a, temp(i)-(day_start-1));       %mean anomaly of venus from J2018
[vol1.ecc_anom_nav, vol1.teta_nav] = kepler1(vol1.manom,vol1.e);
vol1.pos_lettura = fix(vol1.teta_nav/riso)+1;
vol1.pnav(1,1) = vol1.pos(vol1.pos_lettura,1);
vol1.pnav(2,1) = vol1.pos(vol1.pos_lettura,2);
vol1.pnav(3,1) = vol1.pos(vol1.pos_lettura,3);

vol1.plot_posN = plot3(vol1.pnav(1,1),vol1.pnav(2,1),vol1.pnav(3,1),'magenta .', 'MarkerSize',5); hold on;