# Comparing DIC and WAIC for Multilevel Models with Missing Data
1. employee.dat provides the real data example in the paper.
2. model1.imp - model3.imp are the Blimp code for computing likelihoods. model1.blimp-out - model3.blimp-out are the Blimp outputs. ll.model1.csv - ll.model3.csv are the save likelihood results
3. In analyze.R, we extract the likehoods in ll.model1.csv - ll.model3.csv to compute DIC1, DIC2, and WAIC.
4. MC error.csv provides the Monte Carlo errors of each index in the simulation.
