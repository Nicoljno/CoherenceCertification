# Coherence certification methods
Coherence certification by semidefinite programming hierarchies and linear programming.

# List of folders

- Dusseldorf hierarchy

- BlockMatrixOptimization
  - dualBMMl1.m: function to compute the dual of the SDP relaxation at level 1;
  - minwitness.m: function to compute the eigenvectors which minimise the witness after passing through a quantum channel;
  - optimization.m: main file for the optimization search.
  - membership_test.m: Block-moment-matrix relaxation for obtaining upper bounds on the critical visibility for equiangular tight frames (ETF) of dimension dim.
  - ETF_2d.jl: Julia file to obtain 2d equiangular tight frames in dimension d
  - ETF_2d.mat: saved equiangular tight frames used in membership_test.m

- LinearProgramOptimization
  - Binarization_JM_LP.jl: Julia file used to estimate the critical thresholds for joint measurability of the binarization of a given set of states
  - results_depol_damp.txt: saved results for composition of phase damping and depolarising channel
  - results_depol_deph.txt: saved results for composition of dephasing and depolarising channel

