# A Distributed-Memory simulation of the Andromeda-Milky Way collision.
## What is this?
The Andromeda and Milky Way galaxies will collide in approximately 5 billion years, forming Milkdromeda or the Andromeda-Milky Way galaxy. 
This is a distributed memory simulation of this collision run on the ''Supercomputer''/HPC system at the University of Melbourne done as part of coursework for Masters. 

## How was this done?
OpenMPI was used for message passing on the supercomputer/HPC system, work was split evenly among the number of worker nodes and each worker node ran 
an implementation of the Barnes-Hut algorithm 

## I wanna read more
Thanks for your interest, I am pretty proud/happy of this implementation, my report is available [here](report.pdf)

## Caveats
Well I do not simulate general relativity, this is pretty simple Newtonian physics, which is mostly fine
apart from the fact there are two giant black holes at the middle of each galaxy. So this simulation is really a very rough approximation 
of what would happen, to get a more accurate solution, simulation of general relativity is required. 

## Neat facts about this implementation
I do not use Runge-Kutta! Yup, I use something called leapfrog integration which keeps a bound on the error. 
This is critical since we need to conserve momentumn (because you know Physics be like that), something Runge-Kutta cannot do.



