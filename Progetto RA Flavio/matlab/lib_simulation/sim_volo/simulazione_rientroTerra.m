%% CODICE SIMULAZIONE ORBITA RIENTRO SOI TERRA

simF = struct;

[simF.p_start1T, simF.ind_start1T, ...
 simF.p_start2T,simF.ind_start2T] = punto_cont(final.orb,final.orbC,26606,26995, 1, 62832);

simF.E = 2*atanh(sqrt((final.e-1)/(final.e+1))*tan(final.Tet(1)*0.5));                                   %anomalia eccentrica iniziale
simF.M0 = final.e*sinh(simF.E)-simF.E;                                                                                 %anomalia media iniziale
if simF.M0 > pi
    simF.M0 = simF.M0 - 2*pi;
end

simF.n = sqrt(terra.M*6.67e-20/(-final.a^3));                                                                        %moto medio
simF.time_sim = time_sim;                                                                                                         %tempo velocità simulazione [s]

%simulazione
simF.M = simF.M0;
simF.tet = -5;
simF.time = 0;

while simF.tet < 0
    [~,simF.tet] = kepler2(simF.M,final.e);                                                                             %calcolo angolo anom_true
    if simF.tet > pi
        simF.tet = simF.tet -2*pi;
    end
    simF.r = final.p/(1+final.e*cos(simF.tet));
    simF.pos = final.Meulero*(simF.r.*[cos(simF.tet); sin(simF.tet); 0]);
    if simF.tet >0
        break;
    end
    simF.plot = plot3(simF.pos(1,1), simF.pos(2,1), simF.pos(3,1), 'magenta .', 'MarkerSize', 8); hold on
    pause(0.1)
    delete(simF.plot)
    simF.time = simF.time + simF.time_sim;
    simF.M = simF.M0 + simF.n*simF.time;
end


%% Accensione motori per finire sull'orbita di parcheggio
simF.tempo2 = simF.time;
simF.tet0 = simF.ind_start2T*riso;
simF.n = sqrt(terra.M*6.67e-20/(final.rp^2));

simF.time_sim = 10;
simF.tet = simF.tet0;
simF.num_giri = 2;

while simF.tet < simF.tet0 + simF.num_giri*2*pi
    simF.pos = final.Meulero*(final.rp.*[cos(simF.tet); sin(simF.tet); 0]);
    simF.plot = plot3(simF.pos(1,1), simF.pos(2,1), simF.pos(3,1), 'magenta .', 'MarkerSize', 8); hold on;
    pause(0.1)
    delete(simF.plot)
    simF.time = simF.time + simF.time_sim;
    simF.tet = simF.tet0 + simF.n*(simF.time-simF.tempo2);
end
