# bindingsims
Little Sims of Binding Curves

This is a little repository of code to simulate and visualize binding curves in JS.

The front end applets are shamelessly taken from [ProtSim](https://pubs.acs.org/doi/10.1021/acsomega.2c00560)
but the backend math in `update.js` and `simulation.js` are my own implementation to test out 
different binding curve models. 

The current version this supports are: 
1. Protein-ligand binding (quadratic form)
2. 2 Competitive Ligands binding to 1 protein (cubic form)
3. 2 Proteins binding to a single ligand (cubic form)
4. Protein homodimerizaiton

Next up to do on this: 
1. Enzyme kinetics
2. Cheng-Prusoff?
3. Implement these in python?