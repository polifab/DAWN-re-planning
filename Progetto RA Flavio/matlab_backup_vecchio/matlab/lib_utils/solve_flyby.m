%script creato per cercare l'orbita per l'ultimo flyby
Vp = [-34.191538383     7.8937751956    2.08150145676]';        %velocità del pianeta nel suo nodo target

V_inf_mod = 4.6756024;              %modulo della velocità di ingresso SOI richiesta!

V_inf_min = [v_prova];
c = sqrt(V_inf_mod^2-v_prova^2);
V_inf_min = [V_inf_min; sign_y*c; 0];                    %creato il vettore velocità rispetto venere nella SOI

Vs_min = V_inf_min + Vp;
C = [venere.Meulero [venere.pos(53260,1); venere.pos(53260,2); venere.pos(53260,3)];...
        [0 0 0]                 1                                                                                                       ];
pos_finaleD = C * [venere.SOI*cos(tet_out); venere.SOI*sin(tet_out); 0; 1];
pos_finale = [pos_finaleD(1,1); pos_finaleD(2,1); pos_finaleD(3,1)];

%ho v0 e r0, posso quindi creare l'orbita
[navi.a,navi.e,navi.i_orb,navi.W,navi.w,navi.t_anom,navi.p,navi.ra,navi.rp] = rv_to_orb(pos_finale, Vs_min, venere.u_sol*1e-9);

 M = eulerT(navi.W,navi.w,navi.i_orb);
    
    navi.pos = zeros(Nl,3);
    for j = 1:Nl
        r = navi.p/(1+navi.e*cos(teta(j)));               %modulo del raggio posizione
    
        X = r*cos(teta(j));
        Y = r*sin(teta(j));
        
        navi.pos(j,1) = M(1,1)*X + M(1,2)*Y;
        navi.pos(j,2) = M(2,1)*X + M(2,2)*Y;
        navi.pos(j,3) = M(3,1)*X + M(3,2)*Y;
    end
    
    %calcolo dell'anomalia media zero nel giorno di uscita dal flyby (6
    %Ottobre 2019, giorno numero 646 J2018)
    navi.T_anom0 = 4.6093; %rad
    navi.ecc_anom0 = 2*atan(sqrt((1-navi.e)/(1+navi.e))*tan(navi.T_anom0/2));       %in rad anomalia eccentrica
    navi.mean_anom0 = navi.ecc_anom0 - navi.e*sin(navi.ecc_anom0);              %anomalia media iniziale
    navi.n = sqrt(venere.u_sol*1e-9/(navi.a^3));
    
    %plot_orbP =
    %plot3(navi.pos(:,1),navi.pos(:,2),navi.pos(:,3),'magenta');hold on;               %usato solo in fase di progetto
    %plot3(pos_finale(1,1),pos_finale(2,1),pos_finale(3,1),'magenta .', 'MarkerSize', 9);       %usato solo in fase di progetto
    
    %P = sqrt(2*pi*pi*a^3/(venere.u_sol*1e-9))          %usato solo in fase di progetto
    %disp(["Periodo orbitale = " num2str(P/(24*3600)) "giorni"])     %usato solo in fase di progetto