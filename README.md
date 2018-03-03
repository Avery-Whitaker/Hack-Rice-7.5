# Optimizing Parallelism Overhead in Tree-shaped Computations Using Dependency Programming Models
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
We attempted three approaches to analyzing the structures
- Static analysis with heuristics
- Static analysis with machine learning
- Runtime analysis of computational overhead and structure


## Top Down vs Bottom Up
The optimizations for top-down and bottom-up are very different. We hope to have a general solution to both problems.


## Tests
We develop several measures of success using demo applications
- Basic Full Binary Tree With Constant Work (Top-down and bottom-up)
- Unbalanced Binary tree with variable workloads (As determined by a work function, both top-down and bottom-up)
- Inner Product Computation in M.A.D.N.E.S.S. computation. (Vectors of tree computations)


## Results
We have developed a comprehensive measure of efficiency for basic approaches to chunking tree computations. We have investigated several approaches to tackling this problem. We are in the process of developing preliminary results comparing different approaches to the chunking.
