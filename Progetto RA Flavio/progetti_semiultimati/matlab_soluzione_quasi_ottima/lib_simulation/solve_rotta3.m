%script creato per cercare l'orbita per l'ultimo flyby
Vp = [-34.191538383     7.8937751956    2.08150145676]';        %velocità del pianeta nel suo nodo target

V_inf_mod = 4.6756024;              %modulo della velocità di ingresso SOI richiesta!

V_inf_p = [v_prova];
c = sqrt(V_inf_mod^2-v_prova^2);
V_inf_p = [V_inf_p; sign_y*c; 0];                    %creato il vettore velocità rispetto venere nella SOI

Vs_p = V_inf_p + Vp;
pos_finale = [venere.pos(53260,1); venere.pos(53260,2); venere.pos(53260,3)];

%ho v0 e r0, posso quindi creare l'orbita
[navi3.a,navi3.e,navi3.i_orb,navi3.W,navi3.w,navi3.t_anom,navi3.p,navi3.ra,navi3.rp] = rv_to_orb(pos_finale, Vs_p, venere.u_sol*1e-9);

 M = eulerT(navi3.W,navi3.w,navi3.i_orb);
    
    navi3.pos = zeros(Nl,3);
    for j = 1:Nl
        r = navi3.p/(1+navi3.e*cos(teta(j)));               %modulo del raggio posizione
    
        X = r*cos(teta(j));
        Y = r*sin(teta(j));
        
        navi3.pos(j,1) = M(1,1)*X + M(1,2)*Y;
        navi3.pos(j,2) = M(2,1)*X + M(2,2)*Y;
        navi3.pos(j,3) = M(3,1)*X + M(3,2)*Y;
    end
    
    %calcolo dell'anomalia media zero dalla data di uscita dal flyby
    navi3.T_anom0 = (16899-1)*ris;                                                                               %rad
    navi3.ecc_anom0 = 2*atan(sqrt((1-navi3.e)/(1+navi3.e))*tan(navi3.T_anom0/2));        %in rad anomalia eccentrica
    navi3.mean_anom0 = navi3.ecc_anom0 - navi3.e*sin(navi3.ecc_anom0);                      %anomalia media iniziale
    navi3.n = sqrt(venere.u_sol*1e-9/(navi3.a^3));
    
    %plot_orbP =
    %plot3(navi3.pos(:,1),navi3.pos(:,2),navi3.pos(:,3),'magenta');hold on;               %usato solo in fase di progetto
    %plot3(pos_finale(1,1),pos_finale(2,1),pos_finale(3,1),'magenta .', 'MarkerSize', 9);       %usato solo in fase di progetto
    
    % = sqrt(4*pi*pi*nav3i.a^3/(venere.u_sol*1e-9))          %usato solo in fase di progetto
    %disp(["Periodo orbitale = " num2str(P/(24*3600)) "giorni"])     %usato solo in fase di progetto