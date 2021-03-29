# System Identification for Iterative Learning Control
Real-time nonlinear parameter identification for an industrial robot. Code includes measurements and algorithm for a single batch optimization.
This contribution is submitted to IROS 2021. The code features three entry points:

**MASTER_identification** \
Starts the identification of parameters given a recorded movement.

**MASTER_create_symbolic_robot** \
Computes the symbolic solution of the robot model, including inertia, gravitational load and coriolis matrix.
Only required if the robot model is changed.


# Requirements
All algorithms run in MATLAB.

**MATLAB Optimization Toolbox** 

**MATLAB Systen Identification Toolbox**


# Contact Information

Jonas Weigand \
Researcher at the Technical University Kaiserslautern, Chair of Machine Tools and Control Systems\
and at the German Research Center for Artificial Intelligence, Kaiserslautern.\
ORCID: 0000-0001-5835-3106 \
jonas.weigand@mv.uni-kl.de

Gajanan Kanagalingam \
Researcher at the Technical University Kaiserslautern, Chair of Machine Tools and Control Systems\
gajanan.kanagalingam@gmail.de

Martin Ruskowski \
Professor at the Technical University Kaiserslautern, Chair of Machine Tools and Control Systems\
and at the German Research Center for Artificial Intelligence, Kaiserslautern.\
ORCID: 0000-0002-6534-9057 \
martin.ruskowski@mv.uni-kl.de

March 2021
