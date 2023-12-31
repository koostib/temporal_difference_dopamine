# Temporal difference learning algorithm in the ventral tegmental area dopamine neurons

The first is an IgorPro procedure file called ‘Stetsenko and Koos Code 1.pxp’ used for the mean-field simulations and for generating the stochastic synaptic inputs for the DA neuron model. Instructions for running the code are included in the procedure file.  

The generated synaptic inputs are saved in text files called ‘T2_input.dat’ , ‘T3_input.dat’ and NMDA_input.dat’ and these are read by XPPAUT. 

 The XPPAUT code ‘ CELLMOD.ode’ simulates the DA cell responses to the  inputs generated by the Igor code and given in the 3 data files.  The data files simulate 2 consecutive 5s trials as described in the header of the code file.  Default parameters in ‘CELLMOD.ode’ are for simulating the omission of a large fully predicted reward. For instructions on running XPPAUT under MacOS please see https://sites.pitt.edu/~phase/bard/bardware/xpp/xpp.html). 

 Paper: "Neuronal implementation of the temporal difference learning algorithm in the midbrain dopaminergic system"

