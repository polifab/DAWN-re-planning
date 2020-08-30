                                Done
EarthMars_orbit.m:
    - computes planets position/velocity, spacecraft velocity,
    v-infty, both at departure and at arrival
    - computes time of flight
    - computes orbital elements of interplanetary trajectory, and some
    of its characteristics
    
orbitAttempt.m:
    %done - tries to compute and visualize cruise orbit around Sun from
    %       initial position and velocity
    - compute position and velocity of the spacecraft from Earth to
    Mars, following the patched conics method
    - time restricted (tf from EarthMars_orbit.m)
    
Earth_park:
    - plots parking orbit around Earth
    - visually good wrt Earth aspect
    - uses 3rd party function 'earth_sphere'
    
planet_sphere:
    - draft
    - ideally the same as earth_sphere but for all planets
    
plot_orbit.m:
    - plots orbit of the planet, given some coes
    
plot_interplanetary.m:
    - plot Earth and Mars orbit
    - plot initial (27/9/07) and final (17/2/09) of Earth and Mars
    - plot spacecraft trajectory for interplanetary travel


                                TODO:
- complete planet_sphere (do later)
%done - plot interplanetary trajectory -> with (bi)ellipse? 
%       time recursive?
- flyby
%done - position of departure from SOI Earth? -----> patched conics: 
%       initial position as planet position
%done - position of arrival to SOI Mars?  -----> patched conics: 
%       initial position as planet positionand 
- departing position from Mars in flyby?
- visualize each substep
- plot together
- extend plot_orbit to all planet with ease
- adapt all files to one simple workflow