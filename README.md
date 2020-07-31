# multiscale-CD4Tcells
About: The files in this repository constitute a multi-approach and multi-scale model of CD4+ T lymphocytes.

Execution: Running the main.m file within the MATLAB environment will start a stochastic simulation of CD4+ T cells responding to the input signals defined in antigen.mat. Note that all the files in the repository must be placed in the same directory as the main.m file during a simulation.

Copyright: Helikar Lab, University of Nebraska–Lincoln.

Licence: GNU General Public License v3.0.

Software requirements:
1. MATLAB R2020a.
2. The Cobra (COnstraint-Based Reconstruction and Analysis) Toolbox is needed to optimise the metabolic fluxes. It is available at https://opencobra.github.io/cobratoolbox/stable/.

Inputs: Within antigen.mat, the user must provide an input signal representing the temporal trajectory of the abundance of pathogens/antigens in each compartment.

Configuration of the main.m file to parameterise the model to describe the problem under consideration:
1. Set up the compartments by providing Qa, Qb, Vt, Vln, and Vc on lines 99 to 103.
2. Specify the initial number of agents in each compartment and the number of CD4+ T cells in the body part corresponding to that compartment: Nt0, Nln0, Nc0, Nt0_all, Nln0_all, and Nc0_all on lines 77 to 82.
3. Specify the initial cytokine concentrations in the compartments: IL20, IL40, IL60, IL120, IL170, IL180, IL210, IL230, IL270, IFNg0, and TGFb0 on lines 119 to 129.
4. Specify the cytokine production rates associated with the inputs: PIL2in, PIL4in, PIL6in, PIL12in, PIL17in, PIL18in, PIL21in, PIL23in, PIL27in, PIFNgin, and PTGFbin on lines 106 to 116.
5. Specify the extent and grain of the simulation by providing the number of time steps (steps on line 48) and the step size (step_size on line 49).
6. Specify the number of times the model must be implemented to achieve robust average results: rounds on line 52. This should be decided by a convergence study.
7. Fix the four calibration parameters in relation to the problem under consideration and any experimental data available for calibration: Ease_restim, Ease_ACAD, Ease_memory, and Ease_relax on lines 89 to 92.
