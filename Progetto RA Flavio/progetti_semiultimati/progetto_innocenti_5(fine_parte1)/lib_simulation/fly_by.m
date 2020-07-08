f = struct;
f.day_m = i;              %salvo il giorno in cui sono entrato nella SOI per il flyby
f.Vpmin = [dati.dati_venere(1,2); dati.dati_venere(2,2); dati.dati_venere(3,2)];        %salvo la velocità di venere all'inizio del volo
u_sol = venere.u_sol*1e-9;                                                                                          %mi rispetto sole
pippo = sqrt(u_sol/vol1.p);
f.Vsmin = vol1.Meulero*[-pippo*sin(vol1.teta_nav); pippo*(vol1.e+cos(vol1.teta_nav)); 0];         %velocità navicella inizio volo rispetto sole
f.Vinf_m = f.Vsmin-f.Vpmin;                                                                                               %velocità navicella rispetto sole differenza
f.Vinf_m = inv(vol1.Ct_ven_sol)*[f.Vinf_m;0];
f.Vinf_m = [f.Vinf_m(1,1); f.Vinf_m(2,1); f.Vinf_m(3,1)];                                                  %velocità rispetto il sistema venere

f.pos_m = vol1.pos_start;                                                                                                            %posizione rispetto sistema venere

f.tet_m_ven = vol1.tet_end;
f.u_ven = 6.67e-20*venere.M;

[f.a,f.e,f.i_orb,f.W,f.w,f.t_anom0,f.p,f.ra,f.rp] = rv_to_orb(f.pos_m,f.Vinf_m,f.u_ven);
if (f.t_anom0 > pi)
    f.t_anom0 = (2*pi-f.t_anom0);
end

%poichè orbita retrograda, parto da -f.t_anom0
f.Tet = -f.t_anom0:riso:f.t_anom0;

f.Meuler = eulerT(f.W,f.w,f.i_orb);
[f.sizem,f.sizec] = size(f.Tet);
for z = 1:f.sizec
    f.r = f.p/(1+f.e*cos(f.Tet(z)));
    f.pos_fly = f.r*cos(f.Tet(z));
    f.pos_fly = [f.pos_fly; f.r*sin(f.Tet(z)); 0];              %posizione navicella sistema di riferimento iperbole
    
    f.pos_ven(z,:) = f.Meuler*f.pos_fly;                              %posizione navicella sistema di riferimento venere
end

f.F2 = figure('Name','FlyBy Venere', 'Units','pixel','Position',[0 0 800 600]);
f.p_orb = plot3(f.pos_ven(:,1),f.pos_ven(:,2),f.pos_ven(:,3),'magenta'); hold on;
f.p_orbV = plot3(0,0,0,'green .','MarkerSize',10);
axis equal;
xlabel('x(km)');
ylabel('y(km)');
zlabel('z(km)');

%codice di plot navicella in volo
if simula
    simulazione_flyby
end

%codice creazione nuova traiettoria
close (f.F2)
solve_2_fly


clear u_sol pippo