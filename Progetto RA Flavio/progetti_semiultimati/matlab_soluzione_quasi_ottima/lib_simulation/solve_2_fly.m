%calcolo velocit� uscita navicella
f.tempo_volo = time_hyp(f);
vol2.temp_start = f.day_m+f.tempo_volo/24;
[~,~,~,~,vol2.datiSoi] = pianeti_pos(J2018+vol2.temp_start,mercurio,venere,terra,marte,riso);
posi_ven = [vol2.datiSoi.dati_venere(1,1); vol2.datiSoi.dati_venere(2,1); vol2.datiSoi.dati_venere(3,1)];
f.Ct_2 = [eye(3,3), posi_ven; 0 0 0 1]*[eulerT(0,dati.tanom_ven,0), zeros(3,1); 0 0 0 1]*...
               [venere.Meulero, zeros(3,1); 0 0 0 1];
vol2.u = 6.67e-20*1.9891e30;           

f.zz = sqrt(f.u_ven/f.p);
f.Vinf_p = f.Meuler*[-f.zz*sin(f.Tet(f.sizec)); f.zz*(f.e+cos(f.Tet(f.sizec))); 0];                                 %riferito rispetto venere
f.Vinf_p = f.Ct_2*[f.Vinf_p; 0];                                                                                                            %riferito rispetto il sole (non finale)
f.Vinf_p = [f.Vinf_p(1,1); f.Vinf_p(2,1); f.Vinf_p(3,1)];                                                                      %velocit� uscita soi rispetto frame sole

f.Vp_p = [vol2.datiSoi.dati_venere(1,2); vol2.datiSoi.dati_venere(2,2); vol2.datiSoi.dati_venere(3,2)];                    %velocit� venere rispetto sole

f.Vs_p = f.Vinf_p+f.Vp_p;
vol2.vel_start = f.Vs_p;                                                                                                                         %velocit� uscita navicella soi rispetto sole

%calcolo posizione uscita navicella
vol2.pos_start_ven = f.pos_ven(f.sizec,:)';                                                                                           %posizione uscita soi rispetto venere
vol2.pos_start = f.Ct_2*[vol2.pos_start_ven; 1];                                                                                 %posizione uscita soi rispetto sole(non finale)
vol2.pos_start = [vol2.pos_start(1,1); vol2.pos_start(2,1); vol2.pos_start(3,1)];                                %posizione uscita soi rispetto sole

[vol2.a,vol2.e,vol2.i_orb,vol2.W,vol2.w,vol2.t_anom0,vol2.p,vol2.ra,vol2.rp] = rv_to_orb(vol2.pos_start,vol2.vel_start,vol2.u);

vol2.Meulero = eulerT(vol2.W,vol2.w,vol2.i_orb);
    
vol2.pos = zeros(Nl,3);
for j = 1:Nl
    vol2.r = vol2.p/(1+vol2.e*cos(teta(j)));               %modulo del raggio posizione
  
    X = vol2.r*cos(teta(j));
    Y = vol2.r*sin(teta(j));
        
    vol2.pos(j,1) = vol2.Meulero(1,1)*X + vol2.Meulero(1,2)*Y;
    vol2.pos(j,2) = vol2.Meulero(2,1)*X + vol2.Meulero(2,2)*Y;
    vol2.pos(j,3) = vol2.Meulero(3,1)*X + vol2.Meulero(3,2)*Y;
end
vol2.plot_orbN = plot3(vol2.pos(:,1),vol2.pos(:,2),vol2.pos(:,3),'magenta');hold on

vol2.period = sqrt(4*pi*pi*vol2.a^3/vol2.u);
disp(["Periodo orbitale = " num2str(vol2.period/(24*3600)) "giorni"])     %usato solo in fase di progetto
abc = 1;