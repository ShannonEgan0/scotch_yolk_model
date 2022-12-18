# Scotch Yolk Model

#### External Libraries required:    
    numpy
    matplotlib
## Introduction
A time domain simulation of acceleration and velocity created for a flywheel payload using both a straight crankshaft and a 
scotch yolk. The simulation demonstrates the perfect sinusoidal response created by a scotch yolk, and the deviation from the 
perfect sinusoid induced by a straight crankshaft attached to a flywheel.
Different properties of the simulated setup can be adjusted to yield the acceleration and velocity response yielded by 
adjusting the mounting angle of the scotch yolk, or adjusting the motor/wheel rotation frequency.

The program uses Tkinter for a convenient interface, and displays the simulated output using matplotlib and Funcanim.

## Useful Functions
    ang_velocity(freq)
    Calculates angular velocity from the input frequency, and returns an adjusted time step value to maintain fluid motion of 
    the animation regardless of simulated frequency.
    Returns angular velocity in degrees/second and time_step in seconds

    piston(angular_velocity, radius, time, crank_length=1, initial_angle=0, yolk_angle=0, scotch=False)
    Angular velocity in degrees/s, radius in degrees, time in seconds, crank_length in m, initial_angle in degrees, yolk_angle
    in degrees, scotch as a bool
    Leaving the scotch parameter as false will simulate a rod connected directed to the rim of a flywheel at the specified radius,
    while setting it to True will simulate movement with a scotch yolk mounted to the flywheel.
    Returns flywheel angle, flywheel cartesion position (list: [x_pos, y_pos]) and piston position
    
    inst_v(data, time)
    Calculates the derivative of an input time series data set. This will return velocity from displacement, and acceleration from 
    velocity in turn. Returns a series of len(data) - 1, as no derivative value will be returned for index 1.
    Returns new array of the deriviative of the original
    
## Usage
Each of the parameter fields are to be filled in to dictate the properties of the simulation, noting that yolk angle will have no bearing 
on the straight connection mode, and the Connecting Rod Length will have no impact on the scotch yolk. Description of parameters below:


<img align="center" src="Readme Photos/Menu.jpg" height="150">
<figcaption align="center">Menu example</figcaption>   


| Field                 | Effect                                                                |
|-----------------------|-----------------------------------------------------------------------|
| Crank Radius          | Radius of the crank connection to the rim of the motor flywheel       |
| Motor Frequency       | Motor rotation frequency, or frequency of rotation of the flywheel    |
| Connecting Rod Length | Length of connection between flywheel connection and shaft connection |
| Yolk Angle            | Angle of rotation of the scotch yolk, 0° being vertical               |


<img align="center" src="Readme Photos/Straight Crank.jpg" height="220">
<figcaption align="center">Straight crank connection example</figcaption>   


<img align="center" src="Readme Photos/Scotch Yolk.jpg" height="220">
<figcaption align="center">Scotch Yolk example</figcaption>

## To Do
- Restructure so that velocity and acceleration values are calculated in their own loops
- Improve plotting efficiency for animation
- Include shaft position in scotch yolk for more apparent indication when compared to straight model

## Contact
Shannon Egan - shan@egan.mobi

Link to Project - https://github.com/ShannonEgan0/scotch_yolk_model