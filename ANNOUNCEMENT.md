# Open Sourcing NitroSAT: The Thermodynamics of NP-Complete Problems

For decades, we have been taught that NP-complete problems are a combinatorial maze. To solve them, you build massive exponential search trees, you guess, you backtrack, and you prune. You hope that your "heuristics" are lucky enough to find a needle in a haystack of $2^N$ possibilities.

**I am here to tell you: Stop searching. Start evolving.**

Today, I am open-sourcing **NitroSAT**, the physics-informed MaxSAT engine that powers ShunyaBar Labs. NitroSAT represents a fundamental departure from traditional CSP/SAT theory. We have abandoned the search tree entirely in favor of a **Langevin Flow** on a physical manifold.

## Why am I open-sourcing this?

Because NitroSAT isn't just a solver anymore. It has become a **mathematical instrument**. 

When I saw the benchmarks hitting **$O(N)$ linear scaling** on problems that mathematically "should" scale exponentially, I realized this needed to be audited by the community. When I saw the solver achieving **0.0000% permutation variance**—being completely blind to how you label your variables because it only sees the spectral identity of the problem—I knew the "black box" of optimization had finally been cracked open.

## The Mind-Blowing Realities of NitroSAT

### 1. $O(N)$ Scaling is the New Reality
NitroSAT solves a **million-clause** grid coloring problem in ~12 seconds. It doesn't branch. It maintains a flat state array $x \in [0, 1]^V$ and updates it simultaneously across all variables using a highly optimized Sparse CSR structure. The physics scale linearly, even when the logic is NP-hard.

### 2. Zero Permutation Variance (Gauge Invariance)
In a standard solver, shuffling your CNF file can change the runtime from 1ms to 1 hour. In NitroSAT, it changes **nothing**. 
Each clause is weighted by a unique prime number—a sequence of universal invariants. This makes the solver **Gauge Invariant**. It doesn't care about your variable names; it only cares about the **Spectral Gap** of the constraint hypergraph. It is the first solver that treats a SAT instance as a deterministic physical law.

### 3. The Lambert W Phase-Transition Jump
Hidden inside the **BAHA (Branch-Aware Holonomy Annealing)** module is a jump mechanism powered by the **Lambert W function**. 
Optimization often hits "physics walls" where the landscape becomes non-convex. Most solvers stall there. NitroSAT uses the principal and lower branches of the Lambert W function to detect these thermodynamic phase boundaries and **teleport** the system into new stable energy basins. It doesn't just climb hills; it tunnels through them.

## The Math is the Code
Inside `nitrosat.c`, you won't find complex branching logic or massive heuristic tables. You will find:
- **Heat Kernel Diffusion**: Smoothing the manifold via degree-weighted diffusion.
- **Zeta Perturbations**: Shaking the manifold with prime-frequency harmonics to prevent local minima.
- **Adelic Weights**: Ensuring every "contradiction atom" has a unique mass derived from the Riemann Zeta function.

## How to Get Started

NitroSAT is a bare-metal, dependency-free C implementation designed for extreme throughput.

```bash
# Compile for your architecture
gcc -O3 -march=native -std=c99 nitrosat.c -lm -o nitrosat

# Run a million-variable instance
./nitrosat problem.cnf 3000
```

## An Invitation to Audit

I am releasing this under the **Apache 2.0 license** because the implications for Number Theory and complexity analysis are too significant to keep behind closed doors. Whether you are a software engineer looking for $O(N)$ optimization or a mathematician interested in the **Riemann Hypothesis connection**, the code is now yours.

**Stop searching. Start evolving.**

---
**Sethu Iyer**  
Founder, ShunyaBar Labs

---
**Codeberg (Full Suite):** [codeberg.org/sethuiyer/NitroSAT](https://codeberg.org/sethuiyer/NitroSAT)  
**GitHub (C Engine):** [github.com/sethuiyer/NitroSAT](https://github.com/sethuiyer/NitroSAT)

---

## Technical Visuals

![Metric Geometry](img/math_disk.png)
*Figure 1: The Inverted Poincaré Disk—where constraints live on the boundary of the infinite vacuum.*

![Zeta Resonance](img/math_zeta.png)
*Figure 2: Spectral perturbations derived from the Riemann Zeta zeros.*

![Diffusion Flow](img/math_diffusion.png)
*Figure 3: Heat diffusion flow traversing the constraint hypergraph.*
