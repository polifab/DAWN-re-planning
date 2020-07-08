vol2.pnav = eye(3,1);
vol2.manom = mean_anom_t(vol2.mean_anom0, venere.u_sol ,vol2.a, temp(i)-(day_start2-1));       %mean anomaly of venus from J2018
[vol2.ecc_anom_nav, vol2.teta_nav] = kepler1(vol2.manom,vol2.e);
vol2.pos_lettura = fix(vol2.teta_nav/riso)+1;
vol2.pnav(1,1) = vol2.pos(vol2.pos_lettura,1);
vol2.pnav(2,1) = vol2.pos(vol2.pos_lettura,2);
vol2.pnav(3,1) = vol2.pos(vol2.pos_lettura,3);

vol2.plot_posN = plot3(vol2.pnav(1,1),vol2.pnav(2,1),vol2.pnav(3,1),'magenta .', 'MarkerSize',5); hold on;