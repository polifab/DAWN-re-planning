function [exit_manom] = mean_anom_t(in_manom,u_sol,a,day)
    %periodo di campionamento 24 ore terrestri
    %il giorno che devo ricevere in ingresso è relativo alla
    %"differenza" dal 1 gennaio 2000
    
    %ricorda u_sol è in m^3/s^2, per passare in km^3/h^2 devo
    %moltiplicare per 12.96e-3
    
    u = u_sol*12.96e-3;        %km^3/h^2
    tmp = (day-1)*24;       %h = 1ora
    n = (u/(a^3))^0.5;      %moto medio
    exit_manom = n*tmp+in_manom;

end