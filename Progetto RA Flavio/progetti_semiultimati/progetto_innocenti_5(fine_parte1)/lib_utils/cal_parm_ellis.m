function [a,e,vp,p,period] = cal_parm_ellis(ra,rp,mu)
    a = 0.5*(ra+rp);
    e = (ra-rp)/(ra+rp);
    p = a*(1-e^2);

    vp = sqrt(2*mu*((ra/rp)/(ra+rp)));
    period = (2*pi)*sqrt(a^3/mu);
end