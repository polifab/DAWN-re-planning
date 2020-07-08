%ricevo in ingresso la f_infinito desiderata e la v_inf modulo desiderato
v_infm = 5.33;
f_inf = 100*pi/180;

sfb.e = -1/cos(f_inf);
sfb.u_ven = 6.67e-20*venere.M;
sfb.a = -0.5*(sfb.u_ven/(0.5*v_infm^2-sfb.u_ven/venere.SOI));
sfb.p = -sfb.a*(sfb.e^2-1);

sfb.f_real = acos((1/sfb.e)*((sfb.p/venere.SOI)-1));
sfb.v_infm = [sqrt(sfb.u_ven/sfb.p)*sin(sfb.f_real);
                 -sqrt(sfb.u_ven/sfb.p)*(sfb.e+cos(sfb.f_real));
                 0];
                 
sfb.teta = sfb.f_real:-ris:-sfb.f_real;

[sfb.dimr, sfb.dimc] = size(sfb.teta);
sfb.pos = zeros(sfb.dimc,3);
for i = 1:sfb.dimc
    sfb.r = sfb.p/(1+sfb.e*cos(sfb.teta(i)));
    sfb.pos(i,1) = sfb.r*cos(sfb.teta(i));
    sfb.pos(i,2) = sfb.r*sin(sfb.teta(i));
    
end
plot3(sfb.pos(:,1),sfb.pos(:,2),sfb.pos(:,3));


%calcolo della posizione finale sapendo quella iniziale, è una parte test
sfb.Xf = venere.SOI*cos(sfb.f_real);
sfb.Yf = venere.SOI*sin(sfb.f_real);
sfb.Xv =