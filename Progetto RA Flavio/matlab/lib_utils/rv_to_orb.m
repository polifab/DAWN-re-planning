%ricevo in ingresso i vettori 3x1 di r e v (espressi in km) e ne estrapolo tutti i parametri
%orbitali; mu lo devo ricevere in ingresso come km^3/s^2

function [a,e,i,W,w,teta,p,ra,rp] = rv_to_orb(r,v,mu)
    r0 = sqrt(dot(r,r));
    v0 = sqrt(dot(v,v));
    
    h_v = cross(r,v);                               %VETTORE MOMENTO ANGOLARE
    e_v = (cross(v,h_v)/mu-r/r0);                   %VETTORE ECCENTRICITA
    
    E = 0.5*v0^2-mu/r0;                             %energia orbitale
    a = -mu/(2*E);                                  %semiasse maggiore
    e = sqrt(dot(e_v,e_v));                         %eccentricità
    h = sqrt(dot(h_v,h_v));
    
    alfa = cross([0 0 1]',h_v);
    n_vers = alfa/(sqrt(dot(alfa,alfa)));           %versore normale
    
    %calcolo parametri orbitali
    i = acos(dot([0 0 1]',h_v)/h);                  %inclinazione orbitale
    
    if (dot(n_vers,[0 1 0]') > 0)                   %linea dei nodi
        W = acos(dot([1 0 0]',n_vers));
    else
        W = 2*3.1416 - acos(dot([1 0 0]',n_vers));
    end

    if (dot(e_v,[0 0 1]') > 0)                      %argomento del periasse
        w = acos(dot(n_vers,e_v)/e);
    else
        w = 2*3.1416 - acos(dot(n_vers,e_v)/e);
    end
    
    if (dot(r,v) > 0)                               %anomalia vera
        teta = acos(dot(e_v,r)/(e*r0));
    else
        teta = 2*3.1416 - acos(dot(e_v,r)/(e*r0));
    end
    
    p = h^2/mu;
    ra = a*(1+e);
    rp = a*(1-e);



end