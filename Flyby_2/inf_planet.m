function [inf] = inf_planet(ID_Planet)
%% This funciotn contains the radius, radius of the sphere of influnce and name of all planets
%ID_Planet -> Code that indicates the planet
%Inf_Planet is the matrix that contains the planet's data, the rows contains
%the radius, the radius of phere influence and name of the planet
%inf is a vector that contains the planet's data request from user in this order
%inf = [(Radius of the planet) (Radius of phere of influence) (Planet's name)]
Inf_Planet = [2439, 112500, "Mercury";
              6051, 616400, "Venus";
              6371, 924600, "Earth";
              3389, 577400, "Mars";
              69911, 48223000, "Jupiter";
              58232, 54432000, "Saturn";
              25362, 5179200, "Uranus";
              24624, 86668000, "Neptune"
              1195, 3130000, "Pluton"];
inf = Inf_Planet(ID_Planet,:);
end

