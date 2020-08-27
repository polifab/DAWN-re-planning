%-------------------- M-file for the project ------------------------------

% atan2d_0_360.m            -> Computes arctan(y/x) in [0; 360]°

% atmosphere.m              -> Computes density for altitudes in [0-1000]km

% bisect.m                  -> Evaluates roots with the bisection method

% coe_from_sv.m             -> Computes the classical orbital elements 
%                              from the state vector

% dcm_from_q.m              -> Calculates the direction cosine matrix from
%                              the quaternion

% dcm_to_euler.m            -> Finds the angles of the classical Euler 
%                              sequence from the direction cosine matrix

% dcm_to_ypr.m              -> Finds the angles of the yaw-pitch-roll 
%                              sequence from the direction cosine matrix

% f_and_g.m                 -> Calculates the Lagrange coefficients

% f_and_g_ta.m              -> Calculates the Lagrange coefficients from 
%                              the true anomaly change since t0

% fDot_and_gDot.m           -> Calculates d/dt of the Lagrange coefficients

% fDot_and_gDot_ta.m        -> Calculates d/dt of the Lagrange coefficients
%                              from the true anomaly change since t0

% gauss.m                   -> Gauss method to calculate state vector of a 
%                              body from angles-only observations at three 
%                              closely-spaced times

% gibbs.m                   -> Gibbs method of orbit determination to 
%                              compute R2 velocity

% ground_track.m            -> Plots the ground track of an earth satellite

% heun.m                    -> Predictor-corrector method to integrate 1st
%                              order differential equations

% integrate_thrust.m        -> Numerically integrate Equation 6.26 during 
%                              the delta-v burn

% interplanetary.m          -> Determines the spacecraft trajectory from
%                              the sphere of influence of two planets

% J0.m                      -> Computes the Julian day number at 0 UT for 
%                              any year between 1900 and 2100

% kepler_E.m                -> Newton's method to solve Kepler's equation  
%                              for the eccentric anomaly

% kepler_H.m                -> Newton's method to solve Kepler's equation  
%                              for the hyperbolic eccentric anomaly

% kepler_U.m                -> Newton's method to solve the universal 
%                              Kepler equation for the universal anomaly

% lambert.m                 -> Solves Lambert's problem

% los.m                     -> Determine whether the Earth is in the los 
%                              between the satellite and the Sun

% LST.m                     -> Calculates the local sidereal time

% lunar_position.m          -> Calculates the geocentric equatorial 
%                              position vector of the Moon given the Julian
%                              day

% month_planet_names.m      -> Returns the name of the month and the planet
%                              corresponding to the IDs

% orbit.m                   -> Compute orbit of the spacecraft

% planet_elements_and_sv.m  -> Calculates the orbital elements and the 
%                              state vector of a planet from the date and 
%                              universal time

% q_from_dcm.m              -> Calculates the quaternion from the direction
%                              cosine matrix

% ra_and_dec_from_r.m       -> Calculates the right ascension and the
%                              declination from the geocentric equatorial 
%                              position vector

% rk1_4.m                   -> Runge-Kutta for 1st order differential eqs

% rkf45.m                   -> Runge-Kutta-Fehlberg for 1st order 
%                              differential equations

% rv_from_observe.m         -> Calculates (r,v) vectors of an object from 
%                              radar observations

% rv_from_r0v0.m            -> Computes (r,v) from (r0,v0) and elapsed time

% rv_from_r0v0_ta.m         -> Computes (r,v) from (r0,v0) and change in 
%                              true anomaly

% rva_relative.m            -> Find the position, velocity and acceleration
%                              of B relative to A in the LVLH frame on A

% solar_position.m          -> Calculates the geocentric equatorial 
%                              position vector of the Sun, given the julian
%                              date

% stumpS.m                  -> Evaluates the Stumpff function S(z)

% stumpC.m                  -> Evaluates the Stumpff function C(z)

% sv_from_coe.m             -> Computes the state vector (r,v) from the 
%                              classical orbital elements

% twobody3d.m               -> Solves 3D two-body problems and plots the 
%                              results