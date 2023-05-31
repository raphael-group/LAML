* Version 0.4b (beta, unstable):
    * Fix bug of EM in computing likelihood of trees with polytomies
    * Bypass the case where cvxpy cannot solve on prolem instance in Mstep
* Version 0.3:
    * Stable version of v0.3p
    * Add checkpoints
* Version 0.3p (parallelized, unstable):
    * Topology search with simulated annealing
    * Paralellize NNI operations
* Version 0.2:
    * Topology Search version 1: Simple NNI approach.
        - Can take polytomy trees and will restrict topology search to the branches introduced by randomly resolving the polytomies first. v0.2 will note the maximum likelihood resolved tree in this search process, then take this as the start tree and explore topologies on the whole tree.
        - Topology Search pipeline: 
            1. Pick a branch according to some strategy. 
                - Strategies: 
                Each of these strategies produces a scoring for all branches, and then each branch attempt is selected according to the strategy. A branch attempt is successful if the likelihood (without recomputing branch lengths) is higher than the current likelihood. If a branch attempt is unsuccessful, then we draw a second branch attempt without replacement. This ends when no more branches are available. 
                    - random: will pick a random branch from the set of allowed branches 
                    - vanilla : scored according to comparing the az-partition tags. the should change strategy takes the max of $(d_{ab}, d_{ac})$
                    - shouldchange: scored according to comparing the az-partition tags. the should change strategy takes the max of $(d_{ab} - d_{bc}, d_{ac} - d_{bc})$ 
            2. Greedily pick an NNI which improves the tree likelihood around the selected branch. 
            3. If likelihood improves, accept this branch.
            4. Exit condition: 
                - exit if the improvement between the last round and the current round is smaller than some user-input convergence epsilon parameter AND the current best tree topology has been seen before AND the likelihood has been the same within some threshold K times (where k is defined according to having >95% confidence in picking user-input proportion of bad branches - default 0.2)
                - exit if the number of NNI iterations has exceed the max allowed
* Version 0.1:            
    * Preliminary version: Scipy and EM solvers to optimize branch lengths, nu, and phi given a tree topology.
