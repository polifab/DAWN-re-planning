                                Done
body_sphere:
    - allows to plot a sphere representing an object
    - surface with image of the body
    - from Mercury to Pluto, plus Vesta, Ceres and the Sun
    
park_orbit.m:
    - plots parking orbit around a body
    - visually good wrt Earth aspect
    - only circular orbits for now
    
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
    - uses only the main body (2-body problem hypothesys)
    - tested with InSight mission dates
    
plot_orbit.m:
    - plots orbit of all the planets, plus Vesta, Ceres and the Sun

entire_mission.m:
    - allows for complete plots of the mission from Earth to Mars
	- to be extended to the rest, if possible
    - missing Mars flyby

                                TODO:
- flyby
- departing position from Mars in flyby?
%done - visualize each substep
- adapt all files to one simple workflow
- maybe erase global mu and substitute it with a local variable?
- could erase plot of main body from intpl_orbit
- double check parameters using Nasa fact sheets
- fit planet images to their inclinations
- help sections