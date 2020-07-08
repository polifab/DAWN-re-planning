% Shiva Iyer (25 November, 2013)
% Calculate the Stumpff functions C(x) and S(x)
function [C,S] = stumpff(Z)
    C = []; S = [];
    for (z = Z)
        if (z > 0)
            sx = sqrt(z);
            C(end+1,1) = (1-cos(sx))/z;
            S(end+1,1) = (sx-sin(sx))/sx^3;
        elseif (z < 0)
            sx = sqrt(-z);
            C(end+1,1) = (1-cosh(sx))/z;
            S(end+1,1) = (sinh(sx)-sx)/sx^3;
        else
            C(end+1,1) = 1/2;
            S(end+1,1) = 1/6;
        end
    end
end