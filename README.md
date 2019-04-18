Sample input
./cloning.out 128 10 -252525 0.001 -0.50 2 10 1 8.50 100


./cloning.out N_max t_cloningtime randomseed timestep(dt) S(biasing variable, 0 for unbiased) v_0(speed) tau(not used) soft(not used) epsilon(episilon>\gamma \times density/2 is homogenous phase, epsilon< \gamma \times density/2 is flocked... gamma is set to 2. See Farrel paper PRL for definitions) Number_of_clones


Movies are stored in Snapshot....XYZ [unbiased] and CloneSnapshot....XYZ [biased]

Make using make


