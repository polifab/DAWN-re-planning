%% plot 3d delle orbite e del Sole

r1 =1.0e+08*[2.1450    3.1434   -0.3555];           %punto di partenza
r2 =1.0e+08*[0.3486   -4.2711   -0.1978];           %punto di arrivo

s_sis = plot3(terra.pos(:,1),terra.pos(:,2),terra.pos(:,3),'black' ...
              ,marte.pos(:,1),marte.pos(:,2),marte.pos(:,3),'red' ...
              ,vesta.pos(:,1),vesta.pos(:,2),vesta.pos(:,3),'green' ...
              ,ceres.pos(:,1),ceres.pos(:,2),ceres.pos(:,3),'blue' ...
              ,satellite.pos(:,1),satellite.pos(:,2),satellite.pos(:,3),'magenta' ...
              ,0,0,0,'yellow .', 'MarkerSize',30);hold on;grid on;

plot3(r1(1),r1(2),r1(3),'black x','MarkerSize',7)
plot3(r2(1),r2(2),r2(3),'black o','MarkerSize',7)
          
xlim([-5e8 5e8])
ylim([-5e8 5e8])
zlim([-5e8 5e8])


%% plot pianeti

%riferiti al J2000
mean_anom_ter0 = -0.0433;                               %mean anomaly Earth [rad]  
mean_anom_mar0 = deg2rad(19.41248);                     %mean anomaly Mars  [rad]  
mean_anom_ves0 = deg2rad(3.408879412123165E+02);        %mean anomaly Vesta [rad]
mean_anom_cer0 = deg2rad(6.069621345711466E+00);        %mean anomaly Ceres [rad]  

% day_i=2826;                 %giorni dal 2000-01-01 fino al 2007-09-27
day_i=4631;                   %giorni dal 2000-01-01 fino al 2012-09-05

vel_sim = 7;                  %velocità della simulazione [giorni]

for i=-1:vel_sim:911
    tmp=day_i+i;
    
    %Earth
    ter_t = eye(3,1);
    manom = mean_anom_t(mean_anom_ter0,terra.u_sol,terra.a,tmp);    %calcolo anomalia media nel tempo
    [~,teta_ter] = kepler1(manom,terra.e);                          %calcolo anomalia vera
    pos_lettura = fix(teta_ter/risoluzione)+1;
    ter_t(1,1) = terra.pos(pos_lettura,1);
    ter_t(2,1) = terra.pos(pos_lettura,2);
    ter_t(3,1) = terra.pos(pos_lettura,3);
    
    %Mars
    mar_t = eye(3,1);
    manom = mean_anom_t(mean_anom_mar0,marte.u_sol,marte.a,tmp);
    [~,teta_mar] = kepler1(manom,marte.e);
    pos_lettura = fix(teta_mar/risoluzione)+1;
    mar_t(1,1) = marte.pos(pos_lettura,1);
    mar_t(2,1) = marte.pos(pos_lettura,2);
    mar_t(3,1) = marte.pos(pos_lettura,3);
    
    %Vesta
    ves_t = eye(3,1);
    manom = mean_anom_t(mean_anom_ves0,vesta.u_sol,vesta.a,tmp);
    [~,teta_ves] = kepler1(manom,vesta.e);
    pos_lettura = fix(teta_ves/risoluzione)+1;
    ves_t(1,1) = vesta.pos(pos_lettura,1);
    ves_t(2,1) = vesta.pos(pos_lettura,2);
    ves_t(3,1) = vesta.pos(pos_lettura,3);
    
    %Ceres
    cer_t = eye(3,1);
    manom = mean_anom_t(mean_anom_cer0,terra.u_sol,ceres.a,tmp);
    [~,teta_cer] = kepler1(manom,ceres.e);
    pos_lettura = fix(teta_cer/risoluzione)+1;
    cer_t(1,1) = ceres.pos(pos_lettura,1);
    cer_t(2,1) = ceres.pos(pos_lettura,2);
    cer_t(3,1) = ceres.pos(pos_lettura,3);
    
    %plot dei pianeti
    p3 = plot3(ter_t(1,1),ter_t(2,1),ter_t(3,1),'black .','MarkerSize',16);
    p4 = plot3(mar_t(1,1),mar_t(2,1),mar_t(3,1),'red .','MarkerSize',19);
    p5 = plot3(ves_t(1,1),ves_t(2,1),ves_t(3,1),'green .','MarkerSize',12);
    p6 = plot3(cer_t(1,1),cer_t(2,1),cer_t(3,1),'blue .','MarkerSize',12);
    
    
    legend('Orbita Terra','Orbita Marte','Orbita Vesta','Orbita Ceres','Orbita Satellite', ...
            'Sole','Punto di partenza','Punto di arrivo','Terra','Marte','Vesta','Ceres');hold on;
    
    [p] = gen_date(tmp);
   
    pause(0.25)
    delete(p3)
    delete(p4)
    delete(p5)
    delete(p6)
end