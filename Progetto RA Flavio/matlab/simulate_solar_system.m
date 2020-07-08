function [mercurio, venere, terra, marte] = simulate_solar_system(mercurio,venere,terra,marte,Nl,teta,ris,simula,plot_legend)
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
    time_sim = 10*60;                               %simulazione temporale impostata a 10 minuti, potrebbe essere troppo veloce per alcuni casi
    
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
    J2018 = fix(anno_g*18/(24*3600))+2;     %primo gennaio 2018, il +2 è la correzione attuata per confronto con le efemeridi
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
    
    
    %parte di codice dove ci metto il calcolo delle traiettorie di volo
    %codice relativo alla prima parte di volo
    day_start = 250;
    solve_1_fly;
    flag_volo = 1;
    day_start2 = 1000;
    day_start3 = 2500;
    buf_posZ = 1e10;
    
    
    for i= 230:1250
        
        
        %calcolo posizione pianeti giorno per giorno e li plotto
        [mer_t,ven_t,ter_t,mar_t,datiSoi] = pianeti_pos(J2018+temp(i),mercurio,venere,terra,marte,riso);
  
        
        if i >= day_start && flag_volo == 1
            if i == day_start
                pause(1)
                solve_exit_T
            end
            plot_navicella1
            control_SOI
        end
        if i>= day_start2 && flag_volo == 2
            plot_navicella2
            control_SOI2
        end
        if i>= day_start3 && flag_volo == 3
            if i == day_start3
                solve_rientro
            end
            plot_navicella3
            if vol4.pnav(3,1) < 0
                flag_volo = 4;
                delete(r.plotOrb);
                delete(vol4.plot_posN);
                r.plotOrb = plot3(r.posB(:,1),r.posB(:,2),r.posB(:,3),'magenta');
            end
        end
        if i>= day_start3 && flag_volo ==4
            plot_navicella4
            control_SOI3
        end
        if flag_volo == 5
            close all
            break
        end
        
        
        

        p1 = plot3(mer_t(1,1),mer_t(2,1),mer_t(3,1),'blue .','MarkerSize',9);
        p2 = plot3(ven_t(1,1),ven_t(2,1),ven_t(3,1),'green .','MarkerSize',13);
        p3 = plot3(ter_t(1,1),ter_t(2,1),ter_t(3,1),'black .','MarkerSize',16);
        p4 = plot3(mar_t(1,1),mar_t(2,1),mar_t(3,1),'red .','MarkerSize',19);
        %scrivo la data
        [pt] = gen_date(temp(i));
        if plot_legend
            if (i < day_start || flag_volo >= 4) || ((i >=day_start2 && i<day_start3) || flag_volo == 3)
                leg = legend('orb Terra','orb Venere','orb Marte','orb Mercurio' ...
                        ,'Sole','Mercurio','Venere','Terra','Marte');hold on;
            end

            if i > day_start && flag_volo <= 2
                leg =  legend('orb Terra','orb Venere','orb Marte','orb Mercurio' ...
                        ,'Sole','orb Navicella','Naviella','Mercurio','Venere','Terra','Marte');hold on;
            end

            if i >= day_start3 && flag_volo == 3
                leg =  legend('orb Terra','orb Venere','orb Marte','orb Mercurio' ...
                        ,'Sole','orb Navicella','Naviella','Mercurio','Venere','Terra','Marte');hold on;
            end
        end
        
        
        pause(0.05);
        delete(p1);delete(p2);delete(p3);delete(p4);
        if i >= day_start && i<day_start2
            delete(vol1.plot_posN);
        end
        if i>= day_start2 && flag_volo == 2
            delete(vol2.plot_posN);
        end
        if i>= day_start3 && flag_volo == 3
            delete(vol4.plot_posN);
        end
        if i>= day_start3 && flag_volo == 4
            delete(vol4.plot_posN);
        end
        if plot_legend
            delete(pt);
        end
    end
end