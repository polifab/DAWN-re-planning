                                Done
body_sphere:
    - allows to plot a sphere representing an object
    - surface with image of the body
    - from Mercury to Pluto, plus Vesta, Ceres and the Sun
    
park_orbit.m:
    - plots parking orbit around a body
    - visually good wrt planets aspect
    - only circular orbits for now
    - now allows for orbit inclination
    
gen_orbit.m:
    - computes planets position/velocity, spacecraft velocity,
      v-infty, both at departure and at arrival
    - computes time of flight
    - computes orbital elements of interplanetary trajectory, and some
      of its characteristics
    - for any couple of body from Mercury to Pluto, plus Vesta,
      Ceres and the Sun
       
intpl_orbit.m:
    - plot the interplanetary orbit between two points
    - from position and velocity at departure
    - tof as argument
    - uses the Sun as main body (2-body problem hypothesys)
    
plot_orbit.m:
    - plots orbit of all the planets, plus Vesta, Ceres and the Sun
    - now with year indication, for a better orbit plot

entire_mission.m:
    - allows for complete plots of the mission from Earth to Mars
    - demo of mission after Mars
    - missing Mars flyby
    
escape_hyp.m:
    - plots escape hyperbola
    - from circular parking orbit
    - SOI exit with elements of the interplanetary trajectory
    - spacecraft must start from rightly inclinated parking orbit

                                TODO:
- fit planet images to their inclinations -> body_sphere
- add info about starting planet in intpl_orbit? (section already printed)                                

- flyby
- maybe erase global mu and substitute it with a local variable?
- double check parameters using Nasa fact sheets
- help sections
- generate info for parking orbits
- allow for retrograde escape orbit?
- simplify functions deleting unused arguments
- figure performances
- adapt hyperbola to all planets