# Optimizing Parellism Overhead in Tree-shaped Computations Using Dependency Programing Models
## Hack-Rice-7.5

## Team:
- Avery Whitaker (322 TA)
- Madison Lewis (322 TA)
- Ryan Han (322 TA)
- Tian Lan

## About
Hack Rice 7.5 Project
Uses the Intel CnC Runtime
https://icnc.github.io/


## Background Reading:
- [CnC Programing Model](https://wiki.rice.edu/confluence/download/attachments/4425835/cnc_encyc_TR.pdf?version=1&modificationDate=1326693865955)
- [Intel CnC Documentation](https://icnc.github.io/api/index.html)
- [Use Case for Binary Structure Comptuations (Scientific Function Approximation)](https://charm.cs.illinois.edu/workshops/charmWorkshop2014/slides/CharmWorkshop2014_keynote_robert.pdf)


## Approaches
We attempted three approaches to analyzing the sturcutres
- Static anaylsis with heaustics
- Static anaylsis with machine learning
- Runtime anaylsis of computational overhead and structure


## Top Down vs Bottem Up
The optimizations for top down and bottem up are very different. We hope to have a general solution to both problems.


## Tests
We develop several mesures of success using demo applications
- Basic Full Binary Tree With Constant Work (Top down and bottem up)
- Unbalanced Bianry tree with variable workloads (As determined by a work function, both top down and bottem up)
- Inner Product Computation in M.A.D.N.E.S.S. computation. (Vectors of tree computations)


## Results
We have developed a comprehensive measure of effeceincy for basic appraches to chunking tree computations. We have investigated several approaches to tackling this problem.
