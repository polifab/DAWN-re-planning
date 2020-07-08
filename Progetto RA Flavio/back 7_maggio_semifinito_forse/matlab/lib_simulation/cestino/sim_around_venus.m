Fv =figure('Name','Volo attorno a Venere', 'Units','pixel','Position',[0 0 800 600]);
mu_venere = 6.67e-20*Venere.M;
plot(0,0,'m .','MarkerSize',25); hold on;
axis ([-80e3 1e4, -4.5e4 4.5e4])
xlabel('x(km)');
ylabel('y(km)');


%generazione grafica degli ellissi
orb1 = zeros(N,3);
ra1 = 75e3;                    %km
rp1 = 6551;                     %km
[a1,e1,vp1,p1,period1] = cal_parm_ellis(ra1,rp1,mu_venere);
for i = 1:N
    r = p1/(1+e1*cos(Tet(i)));
    
    orb1(i,1) = r*cos(Tet(i));
    orb1(i,2) = r*sin(Tet(i));
end

orb2 = zeros(N,3);
ra2 = 25e3;                    %km
rp2 = 6551;                     %km
[a2,e2,vp2,p2,period2] = cal_parm_ellis(ra2,rp2,mu_venere);
for i = 1:N
    r = p2/(1+e2*cos(Tet(i)));
    
    orb2(i,1) = r*cos(Tet(i));
    orb2(i,2) = r*sin(Tet(i));
end

orb3 = zeros(N,3);
ra3 = 6.551e3;                    %km
rp3 = 6551;                     %km
[a3,e3,vp3,p3,period3] = cal_parm_ellis(ra3,rp3,mu_venere);
for i = 1:N
    r = p3/(1+e3*cos(Tet(i)));
    
    orb3(i,1) = r*cos(Tet(i));
    orb3(i,2) = r*sin(Tet(i));
end

plot(orb1(:,1),orb1(:,2),'blue',orb2(:,1),orb2(:,2),'red',orb3(:,1),orb3(:,2),'black');hold on;

temp = 0:1:60*24*60;
fly = 1;
fly_c = 0;
for i = 1:10
    for j = 1:60*24*60
        if fly == 1
            if temp(j) < (period1/60)
                n = sqrt(mu_venere/(a1^3));
                [E, Te] = kepler1 (n*60*temp(j), e1);
                pos_lettura = fix(Te/risoluzione)+1;
                p = plot(orb1(pos_lettura,1),orb1(pos_lettura,2),'blue .','MarkerSize',10);
                pause(0.05);
                delete(p);
            else
                fly = 2;
            end
        end
        if fly == 2
            if (temp(j)-(period1/60)) < (period2/60)
                n = sqrt(mu_venere/(a2^3));
                [E, Te] = kepler1 (n*60*(temp(j)-period1/60), e2);
                pos_lettura = fix(Te/risoluzione)+1;
                p = plot(orb2(pos_lettura,1),orb2(pos_lettura,2),'blue .','MarkerSize',10);
                pause(0.05);
                delete(p);
            else
                fly = 3;                
            end
        end
        if fly == 3
            n = sqrt(mu_venere/(a3^3));
            [E, Te] = kepler1 (n*60*(temp(j)-period1/60-period2/60), e3);
            pos_lettura = fix(Te/risoluzione)+1;
            p = plot(orb3(pos_lettura,1),orb3(pos_lettura,2),'blue .','MarkerSize',10);
            pause(0.05);
            delete(p);
        end
    end
end
