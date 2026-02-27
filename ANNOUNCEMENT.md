# Open Sourcing NitroSAT: The Physics-Informed MaxSAT Engine

Imagine submitting a constraint problem so large, so tangled with dependencies and resource conflicts, that traditional solvers simply give up.

You have likely been there. You run a shift scheduling optimization for 1,000 nurses across 15 departments. You configure a microservices dependency validation for a Kubernetes cluster. Or perhaps you are trying to find a Ramsey number for mathematical research. You hit “solve,” and you watch the progress bar stall. The CPU spikes, memory usage climbs, and minutes turn into hours. Eventually, the process times out, crashes, or returns a solution that violates your most critical constraints.

Now, imagine a different reality.

Imagine submitting that same massive, messy problem and getting a near-perfect assignment back in milliseconds. Imagine the solver not only providing the result but also telling you exactly why it couldn’t reach 100% — whether it’s a structural parity contradiction in your logic or a genuine resource conflict that requires human intervention.

That is the reality that **NitroSAT** delivers.

Today, ShunyaBar Labs is proud to announce that we are open-sourcing NitroSAT, the core physics-informed MaxSAT solver that powers the intelligence behind Navokoj’s Pro engine. This isn’t just a code release; it is a thinking shift in how we approach NP-hard optimization. We are moving away from brittle heuristics and search-tree pruning and toward a future where constraint satisfaction is treated as a physical dynamical system.

This release represents resilience, speed, and transparency in an industry often defined by opacity and failure.

## From “Search” to “Physics”

For decades, the dominant approach to solving SAT (Boolean Satisfiability) and MaxSAT problems has been combinatorial search. Solvers like CDCL (Conflict-Driven Clause Learning) operate by navigating a massive binary search tree, encountering conflicts, and learning new clauses to avoid those conflicts in the future.

This approach works well for many problems, but it has a fatal flaw: **Structure**.

When problems are “structured” — when they exhibit deep symmetries, tight variable couplings, or phase-transition boundaries — search trees explode. NitroSAT evolves a solution by treating the constraint hypergraph as a physical system governed by energy, heat, and topology.

### Core Modules

1.  **Heat Kernel Gradient Smoothing**: Applying degree-weighted diffusion operators to smooth the energy landscape, preventing the solver from getting stuck in shallow local valleys.
2.  **Persistent Homology**: Utilizing Betti numbers ($\beta_1$) to detect “holes” or loops in the problem structure (e.g., parity chains in XOR-SAT) and breaking them explicitly.
3.  **Zeta-Guided Resonance Injection**: Injecting “prime harmonics” and golden-ratio perturbations to resonate with the underlying structure, helping the solver tunnel through energy barriers.
4.  **Branch-Aware Holonomy Annealing (BAHA)**: Using the Lambert-W function to detect thermodynamic phase boundaries and calculating “jumps” to new energy basins.

## The Zero-Tuning Approach

NitroSAT relies on fundamental mathematical constants (primes, the golden ratio, eigenvalues of the graph Laplacian), allowing the same parameters to work across wildly different domains without manual tweaking.

### Performance Highlights

-   **N-Queens (25x25)**: 100% in 3.7 seconds.
-   **Exact Cover**: 100% in 4ms.
-   **Planted 3-SAT (α=4.26)**: 100% in 5ms.
-   **Hamiltonian Cycle**: 99.99% (missing 1 clause).

## UNSAT Awareness: The Mirage Trap

By monitoring “heat capacity” (energy function variance) and spectral properties, NitroSAT can detect when a solution is physically impossible (a “Mirage Trap”). This transforms the solver from a simple optimizer into a diagnostic engine that quantifies structural impossibility.

## Global Benchmarks

Across 358 benchmark instances:
-   **Average Clause Satisfaction**: 99.58%
-   **Instances Solved at 100%**: 115 (32.1%)
-   **Solved ≥ 99%**: 340 (95.0%)
-   **Speed**: Fastest 10K+ clause solve in 33ms.

## How to Get Started

NitroSAT is lightweight and dependency-free.

### Compile (C)
```bash
gcc -O3 -march=native -std=c99 nitrosat.c -lm -o nitrosat
./nitrosat <cnf-file> [max-steps]
```

### Script (Lua)
```lua
local nitrosat = require("nitrosat")
local solver = nitrosat.NitroSat.new(instance, { seed = 42 })
local success, steps, sat_count = solver:solve()
```

## Why Open Source?

Optimization is filled with black-box solvers. By open-sourcing NitroSAT under the Apache 2.0 license, we invite the community to audit the math, reproduce the benchmarks, and extend the framework.

The future of optimization is not searching; it is physics.

---
**Codeberg (Full Suite):** [codeberg.org/sethuiyer/NitroSAT](https://codeberg.org/sethuiyer/NitroSAT)  
**GitHub (C Engine):** [github.com/sethuiyer/NitroSAT](https://github.com/sethuiyer/NitroSAT)

---

## Technical Visuals (Preview)

![Adelic Manifold](announcement_1.png)
*Figure 1: Visualization of the Adelic Manifold and Spectral Heat Kernels.*

![Phase Transition](announcement_2.png)
*Figure 2: Real-time detection of thermodynamic phase transitions (The Mirage Trap).*

![Heat Kernel Diffusion](announcement_3.png)
*Figure 3: Multi-scale relaxation via heat kernel diffusion operators.*
