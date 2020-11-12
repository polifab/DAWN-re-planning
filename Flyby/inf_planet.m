function [inf] = inf_planet(ID_Planet)
% INF_PLANET contains radius, radius of the sphere of influence
%   and name of all planets.
%
%   ID_Planet  - code that indicates the planet
%
%   Inf_Planet - matrix containing the planet data
%
%   inf        - contains the planet data request from the user:
% inf = [(Radius of the planet) (Radius of sphere of influence) (Planet name)]
%

    Inf_Planet = [2439, 112500, 22032, "Mercury";
                  6051, 616400, 324858, "Venus";
                  6371, 924600, 398600, "Earth";
                  3389, 577400, 42828, "Mars";
                  69911, 48223000, 126711995, "Jupiter";
                  58232, 54432000, 37939519, "Saturn";
                  25362, 5179200, 5780158, "Uranus";
                  24624, 86668000, 6871307, "Neptune"
                  1195, 3130000, 830 "Pluton"];
    inf = Inf_Planet(ID_Planet,:);

end

