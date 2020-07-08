function [mercurio, venere, terra, marte] = simulate_solar_system(mercurio,venere,terra,marte,Nl,teta,ris)
    %definizione costanti
    L_mercurio = 252.25*(pi/180);       %radianti
    w_mercurio = 77.46*(pi/180);       %argomento del periasse
    
    L_venere = 181.97973*(pi/180);       %radianti
    w_venere = 131.53298*(pi/180);       %argomento del periasse
    
    L_terra = 100.46435*(pi/180);       %radianti
    w_terra = 102.94719*(pi/180);       %argomento del periasse
    
    L_marte = 355.45332*(pi/180);       %radianti
    w_marte = 336.04084*(pi/180);       %argomento del periasse
    
    riso = ris;
    
    %calcolo anomalie medie dei pianeti in J2000
    mean_anom_mer0 = L_mercurio-w_mercurio;   %mean anomaly of mercury J2000
    [ecc_anom_mer0, teta_mer0] = kepler1(mean_anom_mer0,mercurio.e);
    mercurio.m_anom0 = mean_anom_mer0;
    
    mean_anom_ven0 = L_venere-w_venere;   %mean anomaly of venus J2000
    [ecc_anom_ven0, teta_ven0] = kepler1(mean_anom_ven0,venere.e);
    venere.m_anom0 = mean_anom_ven0;
    
    mean_anom_ter0 = L_terra-w_terra;   %mean anomaly of earth J2000
    [ecc_anom_ter0, teta_ter0] = kepler1(mean_anom_ter0,terra.e);
    terra.m_anom0 = mean_anom_ter0;
    
    mean_anom_mar0 = L_marte-w_marte;   %mean anomaly of mars J2000
    [ecc_anom_mar0, teta_mar0] = kepler1(mean_anom_mar0,marte.e);
    marte.m_anom0 = mean_anom_mar0;
    

    %ora che ho le posizioni iniziali in J2000, mi dovrò riferire alla data
    %1 Gennaio 2018 e poi parto col plot
    anno_g = 31557600;                      %anno giuliano in secondi
    J2018 = fix(anno_g*18/(24*3600))+1;     %primo gennaio 2018
    %temporale con campionamento giornaliero
    temp = 0:1:2000;
    %Fs =figure('Name','Sistema solare', 'Units','pixel','Position',[0 0 1024 768])
    s_sis = plot3(terra.pos(:,1),terra.pos(:,2),terra.pos(:,3),'black' ...
                  ,venere.pos(:,1),venere.pos(:,2),venere.pos(:,3),'green' ...
                  ,marte.pos(:,1),marte.pos(:,2),marte.pos(:,3),'red' ...
                  ,mercurio.pos(:,1),mercurio.pos(:,2),mercurio.pos(:,3),'blue' ...
                  ,0,0,0,'yellow .', 'MarkerSize',30);hold on;

    axis equal;
    xlabel('x(km)');
    ylabel('y(km)');
    zlabel('z(km)');
    
    %vettore usato per studiare le velocità della navicella all'ingresso di
    %venere
    v_orbitN = zeros(400,8);
    
    %codice relativo alla ricerca dei flyby, mi da l'ultima rotta
    navi = struct;
    on = 1;
    if on
        tet_out = 0;
        v_prova = 1.351;
        sign_y = 1;
        solve_flyby
        on = 0;
    end
    
    for i= 720:900
        
        %calcolo posizione pianeti giorno per giorno e li plotto
        [mer_t,ven_t,ter_t,mar_t] = pianeti_pos(J2018+temp(i),mercurio,venere,terra,marte,riso);

        p1 = plot3(mer_t(1,1),mer_t(2,1),mer_t(3,1),'blue .','MarkerSize',9);
        p2 = plot3(ven_t(1,1),ven_t(2,1),ven_t(3,1),'green .','MarkerSize',13);
        p3 = plot3(ter_t(1,1),ter_t(2,1),ter_t(3,1),'black .','MarkerSize',16);
        p4 = plot3(mar_t(1,1),mar_t(2,1),mar_t(3,1),'red .','MarkerSize',19);
        %scrivo la data
        [pt] = gen_date(temp(i));
        %pt = title(temp(i));
        legend('orb Terra','orb Venere','orb Marte','orb Mercurio' ...
                    ,'Sole','Mercurio','Venere','Terra','Marte','Satellite');hold on;
        
        
        
        % codice posizione spaziale navicella
        %solve_1_fly
        
        if i >= 659
            if on == 0
                plot_orbN = plot3(navi.pos(:,1),navi.pos(:,2),navi.pos(:,3),'magenta');hold on;
                on = 1;
            end
            
            manom = mean_anom_t(navi.mean_anom0,venere.u_sol,navi.a,i-659+1);       %mean anomaly of venus from J2018
            [ecc_anom_ven, teta_nav] = kepler1(manom,navi.e);
            pos_lettura = fix(teta_nav/ris)+1;
            plot_N = plot3(navi.pos(pos_lettura,1),navi.pos(pos_lettura,2),navi.pos(pos_lettura,3),'magenta .', 'MarkerSize', 10);hold on;
            pause(0.05);
            delete(plot_N);
        end
        
        
        pause(0.01);
        delete(p1);delete(p2);delete(p3);delete(p4);
       % delete(plot_orbN);
        delete(pt);
    
    end
end