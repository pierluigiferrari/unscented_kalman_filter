# Unscented Kalman Filter

This is a simple 2D Unscented Kalman Filter (UKF) implementation in C++. In addition to the implementation of the UKF itself, which is contained in ukf.cpp and ukf.h, main.cpp contains a small program for a visual demonstration of the filter using a simulator that is linked below.

The UKF is used to estimate the state of a moving object using noisy LIDAR and RADAR measurements. For an introduction to UKFs, please refer to the paper ["The Unscented Kalman Filter for Nonlinear Estimation"](https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf).

The simulator provides the ground truth state of the object to be tracked and displays the root mean squared error (RMSE) between the filter estimation and the ground truth.

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets): Handles the communication between the program and the simulator.
  * Run either `./install-mac.sh` or `./install-ubuntu.sh`. It is recommended to use one of these scripts to install uWS. Perform a manual installation as explained below only if these scripts don't work for you.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x.
* The simulator which you can download from [here](https://github.com/udacity/self-driving-car-sim/releases).

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF`

Once you launched the executable, simply run the simulator app and select the UKF simulation.
