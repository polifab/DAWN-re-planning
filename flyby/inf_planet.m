function [inf] = inf_planet(ID_Planet)
%% This funciotn contains the radius, radius of the sphere of influnce and name of all planets
%ID_Planet -> Code that indicates the planet
%Inf_Planet is the matrix that contains the planet's data, the rows contains
%the radius, the radius of phere influence and name of the planet
%inf is a vector that contains the planet's data request from user in this order
%inf = [(Radius of the planet) (Radius of phere of influence) (Planet's name)]
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

