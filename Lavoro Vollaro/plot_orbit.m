function plot_orbit(planet_id)
%   planet_id - planet identifier:
%                1 = Mercury
%                2 = Venus
%                3 = Earth
%                4 = Mars
%                5 = Jupiter
%                7 = Uranus
%                8 = Neptune
%                9 = Pluto

    year = [88,
            225,
            365,
            687,
            4333,
            10759,
            30687,
            60190];

    [~, r0, v0, ~] = planet_elements_and_sv(planet_id,2007,1,1,0,0,0);

    pos = [r0];
    for g = 1:year(planet_id)
        %planet position day by day
        [r, ~] = rv_from_r0v0(r0, v0, g*60*60*24);
        pos = cat(1,pos,r);
    end

    plot3(pos(:,1),pos(:,2),pos(:,3),'--')
    
end