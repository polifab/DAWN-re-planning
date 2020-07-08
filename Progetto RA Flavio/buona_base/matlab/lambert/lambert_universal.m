% Shiva Iyer (28 November, 2013)
% Solve Lambert's problem using universal variables. <r1>,<r2> are position
% vectors around the body <mu> at times <t1>,<t2>. If <path> == 1,solve for
% the long path (delta_nu >= pi), else for the short path (delta_nu < pi).
% Return initial and final velocity vectors <v1>,<v2> and iterations <iter>
function [v1,v2,iter] = lambert_universal(mu, t1, r1, t2, r2, path)
    nr1 = norm(r1);
    nr2 = norm(r2);
    mutt = sqrt(mu)*(t2-t1);

    dnu = acos(dot(r1, r2)/(nr1*nr2));
    A = sqrt(nr1*nr2*(1+cos(dnu)));
    if (path == 1)
        A = -A;
    end

    z = 0;
    v1 = zeros(3,1);
    v2 = zeros(3,1);
    for (iter = 1:5000)
        [C,S] = stumpff(z);
        y = abs(nr1+nr2-A*(1-z*S)/sqrt(C));
        x = sqrt(y/C);
        t = x^3*S+A*sqrt(y);
        if (abs(t-mutt) < 1E-6)
            f = 1-y/nr1;
            g = A*sqrt(y/mu);
            gd = 1-y/nr2;
            v1 = (r2-f*r1)/g;
            v2 = (gd*r2-r1)/g;
            return;
        end

        if (abs(z) > 1E-6)
            Cp = (1-z*S-2*C)/(2*z);
            Sp = (C-3*S)/(2*z);
            tp = x^3*(Sp-1.5*S*Cp/C)+0.125*A*(3*S*sqrt(y)/C+A/x);
        else
            tp = (sqrt(2)/40)*y^1.5+0.125*A*(sqrt(y)+A*sqrt(0.5/y));
        end

        z = z-(t-mutt)/tp;
    end
end