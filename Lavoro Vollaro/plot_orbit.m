function plot_orbit(obj_id, annus)
%   planet_id - planet identifier:
%                1 = Mercury
%                2 = Venus
%                3 = Earth
%                4 = Mars
%                5 = Jupiter
%                6 = Saturn
%                7 = Uranus
%                8 = Neptune
%                9 = Pluto
%               10 = Vesta
%               11 = Ceres
%               12 = Sun
%
%   annus    - year considered, to allow for a more precise plot

    year = [88 
            225 
            365 
            687 
            4331 
            10747 
            30589 
            59800 
            90560 
            1330 
            1682 
            25];
        
    colors = ["g"          %green
              "m"          %magenta
              "b"          %blue
              "r"          %red
              "#A2142F"    %darker red
              "#7E2F8E"    %purple
              "#4DBEEE"    %darker cyan
              "c"          %(bright) cyan
              "#D95319"    %orange
              "#77AC30"    %darker green
              "#EDB120"    %ochre
              "#D95319"];  %orange, not visible due to Sun orbit dimensions

    [~, r0, v0, ~] = planet_elements_and_sv(obj_id,annus,1,1,0,0,0);

    pos = [r0];
    for g = 1:year(obj_id)
        %planet position day by day
        [r, ~] = rv_from_r0v0(r0, v0, g*60*60*24);
        pos = cat(1,pos,r);
    end

    plot3(pos(:,1),pos(:,2),pos(:,3),'--', 'Color', colors(obj_id))
    
end