# Spacecraft Simulator

This repository is a realistic physics engine for all types of rocket vehicles launching from earth's surface. The intent of building and maintaining this project is to demonstrate/strengthen my computational abilities and to provide an in-house solution for simulation and modeling of custom rockets. Physical accuracy is of utmost importance. To support flight software development, the goal is to also provide an interface for avionics to communicate with the simulation in real-time to provide realistic analysis of behavior in-flight.

## Table of Contents
- [Code Structure](#code-structure)
- [Setting Up the Environment](#setting-up-the-environment)
- [Active Features List](#active-features-list)
- [About Me](#about-me)

## Code Structure

The main program and all source get compiled with the `g++-13` C++ compiler. All header files are placed at the `include/` directory and all source files are included in the `source/` directory. Unit tests are written with `googletest` (included as a submodule) and are placed in the `tests/` folder.

To remain readable and easy to maintain, I'm using a combination of functional and object-oriented programming patterns. Most physics models are stateless and are better suited to be represented by collections of functions. The main component used in the simulation is a state vector containing position, velocity, angular velocity, and attitude quaternion (13 floating points). For data structures explicitly used by physical models, I always use a `struct` to encapsulate the data. For subsystems that require their own interfaces for more complex interactions with the environment, _then_ I will prefer to use a `class` or inheritence pattern.

The philosophy for this project is to use the best-suited programming pattern and does not restrict itself to any singular form.

## Setting Up the Environment

I develop on a Macbook Pro using homebrew for package management. The compiler expects the user to have installed `g++-13`. The `configure_environment.sh` script is designed to build the environment you'll need to replicate my processes and use the tools natively.

## Active Features List

| Name | Capabilities | Status |
| ---- | ------------ | ------ |
| WGS84 Model | Conversion between LLA and ECEF coordinate systems | Implemented |
| Rotating reference frame (earth) | Adds Coriolis and Centrifugal non-inertial forces | Implemented |
| Body-frame attitude tracking | Use body-frame angular rates to integrate attitude over time with a quaternion | Implemented |
| Aerodynamic coefficients tables | Interpolate multi-dimensional datasets for aerodynamic coefficients as functions of speed and angle of attack | To-Do |
| Aerodynamic modeling | Implement aerodynamic forces over body surfaces including lift, drag, and grid fins | To-Do |
| Thrust modeling | Interpolate multi-dimensional thrust tables for custom engine types | To-Do |
| Custom thrust vector mounts | Allow user to define how actuator inputs map to unit-thrust-vectors | To-Do |

## About Me

I'm a software engineer with a degree in physics and a burning passion for anything that can fly. I'm experienced in embedded systems, hardware-in-the-loop testing, and numerical modeling. I have ~3 years of experience modeling complex physical systems and am well-versed in Python and C/C++ programming languages.

The purpose of this project is to demonstrate my understanding of rocket dynamics, earth-based navigation, and numerical modeling. This project goes hand-in-hand with a model rocket program that I'm currently working on in parallel. When `SpacecraftSimulator` reaches a deployable state, it will be used to HIL test the flight systems for my rocket.
