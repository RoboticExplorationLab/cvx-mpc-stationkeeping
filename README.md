# Contingency-Aware Station-Keeping Control of Halo Orbits 

A brief, one-sentence description of what your project does and its primary purpose.

## Table of Contents
- [About](#about)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## About
We present an algorithm to perform fuel-optimal stationkeeping for spacecraft in unstable halo orbits with additional constraints to ensure safety in the event of a control failure. To enhance safety, we enforce a half-space constraint on the spacecraft trajectory. This constraint biases the trajectory toward the unstable invariant manifold that escapes from the orbit away from the planetary body, reducing the risk of collision. We formulate a convex trajectory-optimization problem to autonomously generate impulsive spacecraft maneuvers to loosely track a halo orbit using a receding-horizon controller. Our solution also provides a safe exit strategy in the event that propulsion is lost at any point in the mission. We validate our algorithm in simulations of the three-body Earth-Moon and Saturn-Enceladus systems, demonstrating both low total delta-v and a safe contingency plan throughout the mission.

## Organization
This repository contains its own Julia 1.10.1 environment specified by the Project.toml and Manifest.toml files. 

The directoreis in the project are the following: 
