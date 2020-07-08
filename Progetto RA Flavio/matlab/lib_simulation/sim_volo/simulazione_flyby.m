%nella struttura f mi servono: f.Tet è un vettore di teta che va da - anom0
%a + anom0       FUNZIONANTE

%calcolo di M0 (mean anomaly all'inizio del flyby)
f.sim.E = 2*atanh(sqrt((f.e-1)/(f.e+1))*tan(f.Tet(1)*0.5));
f.sim.M0 = f.e*sinh(f.sim.E)-f.sim.E;
if f.sim.M0 > pi                                                                            %correzione algoritmica
    f.sim.M0 = f.sim.M0 -2*pi;
end

f.sim.n = sqrt(6.67e-20*venere.M/(-f.a^3));                                 %moto medio
f.sim.time_sim = time_sim;                                                               %campionamento della simulazione, settata ora a mezz'ora

%inizio simulazione plot
f.sim.M = f.sim.M0;
f.sim.tet = -5;                                                                                 %inizializzazione algoritmica
f.sim.i = 1;

while f.sim.tet <= f.t_anom0
    [~,f.sim.tet] = kepler2(f.sim.M,f.e);
    if f.sim.tet > pi                                                                        %correzione algoritmica
        f.sim.tet = f.sim.tet-2*pi;
    end
    
    if f.sim.tet > f.t_anom0
        break
    end
    f.sim.r = f.p/(1+f.e*cos(f.sim.tet));
    f.sim.pos_fly = f.sim.r*cos(f.sim.tet);
    f.sim.pos_fly = [f.sim.pos_fly; f.sim.r*sin(f.sim.tet); 0];              %posizione navicella sistema di riferimento iperbole
    
    f.sim.pos = f.Meuler*f.sim.pos_fly;
    f.sim.plot = plot3(f.sim.pos(1,1),f.sim.pos(2,1),f.sim.pos(3,1),'magenta .','MarkerSize',8);hold on;
    pause(0.1);
    delete(f.sim.plot)
    f.sim.M = f.sim.M0+(f.sim.i*f.sim.time_sim)*f.sim.n;
    f.sim.i = f.sim.i+1;
end