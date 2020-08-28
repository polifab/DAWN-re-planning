                                Done
Earth_Mars.m:
    - computes planets position/velocity, spacecraft velocity,
    v-infty, both at departure and at arrival
    - computes time of flight
    - computes orbital elements of interplanetary trajectory, and some
    of its characteristics
    
orbitAttempt.m:
    - tries to compute and visualize cruise orbit around Sun from
    initial position and velocity
    - time restricted
    
Earth_park:
    - plots parking orbit around Earth
    - visually good wrt Earth aspect
    - uses 3rd party function 'earth_sphere'
    
planet_sphere:
    - draft
    - ideally the same as earth_sphere but for all planets
    
                                TODO:
- complete planet_sphere (do later)
- plot interplanetary trajectory -> with (bi)ellipse? time recursive?
- flyby
- position of departure from SOI Earth?
- position of arrival to SOI Mars? and departing from it?
- visualize each substep
- plot together