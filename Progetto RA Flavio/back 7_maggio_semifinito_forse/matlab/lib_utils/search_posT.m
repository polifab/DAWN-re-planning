function [index] = search_posT(teta, array)

[rig, col] = size(array);
contr = teta-pi;
if contr < 0
    contr = contr + 2*pi;
end

for i = 1:rig
    angolo = atan2(array(i,2),array(i,1));    
    
    if angolo < 0
        angolo = angolo + 2*pi;
    end
    
    if (abs(contr-angolo) <= 1e-4)
        index = i;
        break;
    end
end
end