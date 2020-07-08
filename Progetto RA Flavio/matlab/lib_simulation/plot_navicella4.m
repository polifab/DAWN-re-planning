vol4.pnav = zeros(3,1);
vol4.m_anom = mean_anom_t(vol4.m_anom0, venere.u_sol, vol4.a, temp(i) - (day_start3-1));
[~, vol4.t_anom] = kepler1(vol4.m_anom, vol4.e);
vol4.pos_lettura = fix(vol4.t_anom/riso)+1;
vol4.pnav(:,1) = r.posB(vol4.pos_lettura,:)';
vol4.plot_posN = plot3(vol4.pnav(1,1), vol4.pnav(2,1), vol4.pnav(3,1), 'magenta .', 'MarkerSize', 5); hold on;