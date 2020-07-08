function [tvolo] = time_hyp(str)
    %str è la struttura puntata per il calcolo del tempo di volo nella SOI
    F0 = acosh((str.e+cos(str.t_anom0))/(1+str.e*cos(str.t_anom0)));
    tp = sqrt(-(str.a^3)/str.u_ven)*(str.e*sinh(F0)-F0);
    Ff = acosh((str.e+cos(-str.t_anom0))/(1+str.e*cos(-str.t_anom0)));
    tf = tp+sqrt(-(str.a^3)/str.u_ven)*(str.e*sinh(Ff)-Ff);
    tvolo = tf/3600;
end