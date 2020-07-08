volF = struct;
volF.tempo_volo = 0;

%% Plot orbita di parcheggio prima dell'uscita dalla soi (tempo simulazione 1 minuto)

%bisogna cercare l'angolo di contatto tra il periasse dell'iperbole e la
%circonferenza, ricerca
[volF.p_start1, volF.ind_start1,...
 volF.p_start2, volF.ind_start2] = punto_cont(r.pos_ip_ven,r.posC_v,1,1,1,62832);

volF.t_anom0 = volF.ind_start2*riso;                                       %potrebbe essere qui l'errore
volF.time_simul = time_sim;                                                            %simulazione impostata al minuto
volF.n = sqrt(6.67e-20*venere.M/(r.rp_ip_venere^3));        %moto medio circolare attorno venere
volF.numero_giri = 4;                                                          %numero dei giri usati per la simulazione
volF.t_anom = volF.t_anom0;

while volF.t_anom <= volF.numero_giri*2*pi+volF.t_anom0
    volF.pos_fly = r.rp_ip_venere*cos(volF.t_anom);
    volF.pos_fly = [volF.pos_fly; r.rp_ip_venere*sin(volF.t_anom); 0];
    
    volF.plot = plot3(volF.pos_fly(1,1),volF.pos_fly(2,1),volF.pos_fly(3,1), 'magenta .','MarkerSize',6);hold on;
    pause(0.05);
    delete(volF.plot)
    volF.tempo_volo = volF.tempo_volo + volF.time_simul;
    volF.t_anom = volF.t_anom0+volF.n*(volF.tempo_volo);
end


%% Plot uscita dal ramo iperbolico, cambio il tempo di simulazione per velocizzare (tempo simulazione = 10 minuti)

volF.m_anom0 = 0;                                                               %sono sicuro, sono sul periasse dell'iperbole
volF.time_simul = time_sim;                                                    %cambio tempo simulazione a 10 minuti
volF.n = sqrt(6.67e-20*venere.M/(-r.a_ip_venere^3));      %moto medio
volF.m_anom = volF.m_anom0;
volF.r = 0;                                                                             %presetto a zero, mi serve per la ciclicità dell'algoritmo
volF.tempo_in = volF.tempo_volo;

while volF.r <= venere.SOI
    [~,volF.tet] = kepler2(volF.m_anom,r.e_ip_venere);
    if volF.tet > pi
        volF.tet = volF.tet - 2*pi;
    end
    
    volF.r = r.p_ip_venere/(1+r.e_ip_venere*cos(volF.tet));
    volF.pos_fly = [volF.r*cos(volF.tet); volF.r*sin(volF.tet); 0];
    volF.pos_fly = r.Meuler_ip_venere*volF.pos_fly;
    volF.plot = plot3(volF.pos_fly(1,1),volF.pos_fly(2,1),volF.pos_fly(3,1), 'magenta .','MarkerSize',6);hold on;
    pause(0.05);
    delete(volF.plot)
    volF.tempo_volo = volF.tempo_volo + volF.time_simul;
    volF.m_anom = volF.m_anom0 + volF.n*(volF.tempo_volo - volF.tempo_in);    
end