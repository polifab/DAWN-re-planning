%% libreria generazione simulazione volo di arrivo a venere, la durata totale la salvo
% poichè mi serve per poi poter generare il codice di rientro sulla terra

%                       DUTATA VOLO = 2742 min = 45.7 ore = circa 2 giorni!

simV = struct;
simV.tempo_inizio = i;                                          %salvo in memoria il giorno dell'arrivo su venere

%algoritmo oneroso, mi servirebbe nella traiettoria finale, la faccio
%calcolare all'inizio
[simV.p_start1T, simV.ind_start1T, ...
 simV.p_start2T,simV.ind_start2T] = punto_cont(vol3.pos_navC,vol3.pos_navT,24440,24930,5458,5746);

%% volo prima orbita, quella iperbolica d'ingresso ORBITA RETROGRADA
simV.E = 2*atanh(sqrt((vol3.e-1)/(vol3.e+1))*tan(vol3.Tet(1)*0.5));
simV.M0 = vol3.e*sinh(simV.E)-simV.E;                                                           %calcolo anomalia media iniziale
if simV.M0 > pi                                                                            %correzione algoritmica
    simV.M0 = simV.M0 -2*pi;
end
simV.n = sqrt(6.67e-20*venere.M/(-vol3.a^3));                           %moto medio
simV.time_sim = time_sim;                                                                %velocità simulazione, settata a 10 minuto

%simulazione
simV.M = simV.M0;
simV.tet = -5;
simV.time = 0;

while simV.tet < 0
    [~,simV.tet] = kepler2(simV.M,vol3.e);
    if simV.tet > pi                                                                        %correzione algoritmica
        simV.tet = simV.tet-2*pi;
    end
    simV.r = vol3.p/(1+vol3.e*cos(simV.tet));
    simV.pos_fly = simV.r*cos(simV.tet);
    simV.pos_fly = [simV.pos_fly; simV.r*sin(simV.tet); 0];         %posizione piano orbitale navicella
    simV.pos = vol3.Meuler*simV.pos_fly;                             %posizione rispetto venere
    simV.plot = plot3(simV.pos(1,1),simV.pos(2,1),simV.pos(3,1),'magenta .','MarkerSize',8);hold on;
    pause(0.1);
    delete(simV.plot)
    simV.time = simV.time+simV.time_sim;
    simV.M = simV.M0+simV.time*simV.n;
    
end

%% Volo 2 Prima accensione motori per rallentare, orbita RETROGRADA
simV.tempo_in1 = simV.time;                     %salvo il tempo di inizio volo 2
simV.tet = pi;
simV.M = pi;                                                %lo so, sono sull'apoasse dell'ellisse
simV.M0 = pi;
simV.n = sqrt(6.67e-20*venere.M/(vol3.aB^3));       %moto medio ellisse

while simV.M < 2*pi
    [~,simV.tet] = kepler2(simV.M,vol3.eB);             %calcolo anomalia vera
    
    simV.r = vol3.pB/(1+vol3.eB*cos(simV.tet));
    simV.pos_fly = simV.r*cos(simV.tet);
    simV.pos_fly = [simV.pos_fly; simV.r*sin(simV.tet);0];
    simV.pos = vol3.MeulerB*simV.pos_fly;
    simV.plot = plot3(simV.pos(1,1),simV.pos(2,1),simV.pos(3,1),'magenta .','MarkerSize',8);hold on;
    pause(0.1);
    delete(simV.plot)
    simV.time = simV.time+simV.time_sim;
    simV.M = simV.M0 + simV.n*(simV.time-simV.tempo_in1);
end


%% Volo 3 ORBITA CIRCOLARE RETROGRADA R=500km (dalla superficie di venere)
simV.tempo_in2 = simV.time;                                                                                                             %memorizzzo il tempo iniziale
simV.pos_in = vol3.pos_startC;
[simV.p_start, simV.ind_start1,~,~] = punto_cont(vol3.pos_navC,simV.pos_in,1,62832,1,1);
simV.tet0 = simV.ind_start1*riso;                                                                                                       %anomalia vera iniziale
simV.n = sqrt(6.67e-20*terra.M/(vol3.aC)^3);
simV.ind_giri = 3;
simV.tet = simV.tet0;
simV.r = norm(simV.pos_in);                                                                                                                           %è orbita circolare
simV.time_sim = time_sim;

while simV.tet < simV.tet0 + 2*simV.ind_giri*pi + (simV.ind_start1T*riso-simV.tet0)
    simV.pos_fly = simV.r*cos(simV.tet);
    simV.pos_fly = [simV.pos_fly; simV.r*sin(simV.tet); 0];
    
    simV.pos = vol3.MeulerC*simV.pos_fly;
    simV.plot = plot3(simV.pos(1,1),simV.pos(2,1),simV.pos(3,1),'magenta .','MarkerSize',8);
    pause(0.1)
    delete(simV.plot)
    simV.time = simV.time + simV.time_sim;
    simV.tet = simV.tet0 + simV.n*(simV.time-simV.tempo_in2);    
end

%% Volo 4 VOLO ORBITA CIRCOLARE POLARE RETROGRADA a 500km dalla superficie di Venere

simV.tet0 = simV.ind_start2T*riso;
simV.n = sqrt(6.67e-20*venere.M/(vol3.RT)^3);                           %moto medio
simV.M = simV.tet0;                                                                      
simV.time_in = simV.time;

while simV.M <= 5*2*pi+simV.tet0
    simV.pos_fly = vol3.RT*cos(simV.M);
    simV.pos_fly = [simV.pos_fly; vol3.RT*sin(simV.M); 0];
    
    simV.pos = vol3.MeulerT*simV.pos_fly;
    simV.plot = plot3(simV.pos(1,1),simV.pos(2,1),simV.pos(3,1),'magenta .','MarkerSize',8);
    pause(0.1)
    delete(simV.plot)
    simV.time = simV.time + simV.time_sim;
    simV.M = simV.tet0 + simV.n*(simV.time-simV.time_in);
end

durata_volo = simV.time/60;
disp("Durata volo = " + num2str(durata_volo))