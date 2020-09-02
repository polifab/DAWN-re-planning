                                Done
gen_orbit.m:
    - computes planets position/velocity, spacecraft velocity,
      v-infty, both at departure and at arrival
    - computes time of flight
    - computes orbital elements of interplanetary trajectory, and some
      of its characteristics
    - for any couple of body from Mercury to Pluto, plus Vesta,
      Ceres and the Sun
    
orbitAttempt.m:
    %done - tries to compute and visualize cruise orbit around Sun from
    %       initial position and velocity
    - compute position and velocity of the spacecraft from Earth to
    Mars, following the patched conics method
    - time restricted (tf from EarthMars_orbit.m)
    
park_orbit.m:
    - plots parking orbit around a body
    - visually good wrt Earth aspect
    - only circular orbits for now
    
body_sphere:
    - allows to plot a sphere representing an object
    - surface with image of the body
    - from Mercury to Pluto, plus Vesta, Ceres and the Sun
    
plot_orbit.m:
    - plots orbit of all the planets, plus Vesta, Ceres and the Sun
    
plot_interplanetary.m:
    - plot Earth and Mars orbit
    - plot initial (27/9/07) and final (17/2/09) of Earth and Mars
    - plot spacecraft trajectory for interplanetary travel
    
intpl_orbit.m:
    - plot the interplanetary orbit between two points
    - from position and velocity at departure
    - tof as argument
    - uses only the main body (2-body problem hypothesys)


                                TODO:
%done in body_sphere - complete planet_sphere (do later)
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
%done - extend plot_orbit to all planet with ease
- adapt all files to one simple workflow
- maybe erase global mu and substitute it with a local variable?