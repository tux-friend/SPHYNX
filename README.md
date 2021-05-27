# SPHYNX
Used code for neutron star merger simulations. Adapted version of SPHYNX - an SPH hydrocode by Rubén M. Cabezón and Domingo García-Senz.

This repository contains the code I used for my master thesis project. The code is adapted to the use of neutron star merger simulations.

Source: https://astro.physik.unibas.ch/en/people/ruben-cabezon/sphynx/
For detailed information, see also: https://www.aanda.org/articles/aa/pdf/2017/10/aa30208-16.pdf

- **Relaxation** contains the code for relaxing a initial 3D model
- **Merger** contains the code for performing a merger simulation of two identical neutron stars 
- **1D_density_profiles** - scripts for creating and plotting 1D density profiles for polytropic and piecewise polytropic EOS
- **plotSPYNXdata** - scripts for plotting models, GW signal, stages of coalescence etc.

When using a tabulated EOS a table e.g. from https://stellarcollapse.org/equationofstate needs to be placed into the **eos** folder.
