# Contingency-Aware Station-Keeping Control of Halo Orbits 

## Table of Contents
- [About](#about)
- [Organization](#organization)

## About
We present an algorithm to perform fuel-optimal stationkeeping for spacecraft in unstable halo orbits with additional constraints to ensure safety in the event of a control failure. To enhance safety, we enforce a half-space constraint on the spacecraft trajectory. This constraint biases the trajectory toward the unstable invariant manifold that escapes from the orbit away from the planetary body, reducing the risk of collision. We formulate a convex trajectory-optimization problem to autonomously generate impulsive spacecraft maneuvers to loosely track a halo orbit using a receding-horizon controller. Our solution also provides a safe exit strategy in the event that propulsion is lost at any point in the mission. We validate our algorithm in simulations of the three-body Earth-Moon and Saturn-Enceladus systems, demonstrating both low total delta-v and a safe contingency plan throughout the mission.

## Organization
This repository contains its own Julia 1.10.1 environment specified by the Project.toml and Manifest.toml files. 

The directories in the project are the following: 
- src: algorithm source code
- examples: contains case studies in the Earth-Moon and the Saturn Enceladus systems
- solution_trajectories: contains the solutions for the two example scenerios that are used to generate the escape trajectories
- escape_trajectories: jupyter notebook to generate the escape trajectories
- refs: reference trajectories used for the examples