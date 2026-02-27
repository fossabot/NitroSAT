## NitroSAT Performance Benchmarks  

![benchmarks](img/benchmarks.png)

### 1️⃣ Large Structured Instances  

| Problem (type) | Variables | Clauses | Satisfaction | Time | Hardware |
|-----------------|-----------|---------|--------------|------|----------|
| **Graph Coloring** | – | 650 000 | **100 %** | 4.6 s | Ryzen 5 5600H |
| **Clique Coloring – `cliquecol_80_10_10`** | 4 760 | 354 890 | **100 %** (5/5 seeds) | **3.5 s** | Laptop |
| **Ramsey R(5,5,5)** | – | 402 752 | **100 %** | 7.5 s | Laptop |
| **Parity (CNFgen)** | – | 485 200 | **100 %** | 2.6 s | – |
| **Counting (CNFgen)** | – | 78 760 | **100 %** | 0.5 s | – |
| **Matching (100‑node)** | – | 400 | **100 %** | 21 ms | – |
| **Van der Waerden** | – | 20 | **100 %** | <1 ms | – |
| **Tiling (8×8 grid)** | – | 580 | **99.1 %** | 154 ms | – |
| **Subset Cardinality** | – | 210 | **95.7 %** | 107 ms | – |

---  

### 2️⃣ Real‑World Scheduling Benchmarks  

| Jobs | Slots | Density | Clauses | Satisfaction | Time |
|------|-------|---------|---------|--------------|------|
| 30 | 5 | 0.30 | 1 095 | **100 %** | 28 ms |
| 50 | 6 | 0.20 | 2 522 | **100 %** | 22 ms |
| 100 | 8 | 0.10 | 7 516 | **100 %** | 26 ms |
| 200 | 10 | 0.05 | 21 150 | **100 %** | 85 ms |
| 500 | 10 | 0.03 | 64 570 | **100 %** | 172 ms |
| **1 000** | 10 | 0.02 | 154 760 | **99.99 %** (UNSAT detection) | 225 s |

---  

### 3️⃣ Random 3‑SAT Phase‑Transition (α ≈ 4.26)  

| Variables (n) | Clauses | Seeds | Avg. Sat. | Min | Max | Std. Dev. |
|---------------|---------|-------|-----------|-----|-----|-----------|
| 300 | 1 278 | 50 | **99.65 %** | 99.37 % | 99.84 % | 0.11 % |
| 500 | 2 130 | 20 | **99.64 %** | 99.44 % | 99.81 % | 0.10 % |
| 1 000 | 4 260 | 10 | **99.65 %** | 99.53 % | 99.72 % | **0.06 %** |

*The variance **shrinks** as the instance size grows, evidencing a scale‑invariant robustness.*  

---  

### 4️⃣ Grid‑Coloring Stress Test (Spectral Scaling)  

| Grid Size | Variables | Clauses | Satisfaction | Time |
|-----------|-----------|---------|--------------|------|
| 10 × 10 | 400 | 1 420 | **100 %** | 0.02 s |
| 20 × 20 | 1 600 | 5 840 | **100 %** | 0.07 s |
| 50 × 50 | 10 000 | 37 100 | **100 %** | 0.45 s |
| 100 × 100 | 40 000 | 149 200 | **100 %** | 2.10 s |
| **1 000 × 1 000** | **4 000 000** | **14 992 000** | **100 %** | **475 s** |

---  

### 5️⃣ XOR‑SAT Stress Test (Parity‑Chain Detection)  

| Instance | Clauses | Result | Time | Detected Loops (β₁) |
|----------|---------|--------|------|----------------------|
| XOR (SAT) | 200 | **100 %** | 3.9 ms | 98 |
| XOR (planted) | 8 000 | **100 %** | 9.5 ms | — |
| XOR (UNSAT) | 2 000 | 95.15 % | 802 ms | 2–26 |

---  

### 6️⃣ UNSAT Awareness – “Mirage” Trap  

| Benchmark | Stopping Satisfaction | Observation |
|-----------|----------------------|-------------|
| Mirage 200 | **85.2 %** | Trap Detected (structural impossibility) |
| Mirage 300 | **95.9 %** | Structural Awareness (phase‑transition signal) |

---  

### 7️⃣ Permutation Invariance (Encoding‑Agnostic)  

| Permutations Tested | Perfect Solves | Std. Dev. |
|---------------------|----------------|-----------|
| 20 random variable/ clause permutations | 20 / 20 (100 %) | **0.0000 %** |

---  

### 8️⃣ Quasigroup / Latin‑Square Completion (New CNF Class)  

| Run Type | Instances (SAT + UNSAT) | Avg. Sat. | Std. Dev. | ≥ 99.9 % | Perfect 100 % |
|----------|--------------------------|-----------|-----------|----------|----------------|
| Single‑seed (8 instances) | 8 (4 SAT + 4 UNSAT) | **99.976 %** | – | 8 / 8 | 1 / 8 |
| Seed sweep (40 runs = 8 × 5 seeds) | 40 | **99.960 %** | 0.036 % | 38 / 40 | – |
| SAT runs (40) | – | **99.985 %** | – | 20 / 20 | – |
| UNSAT runs (40) | – | **99.967 %** (plateau) | – | 0 / 20 | – |

---  

### 9️⃣ Global Benchmark Summary (All 358 Instances)  

| Metric | Value |
|--------|-------|
| Total instances evaluated | **358** |
| Solved at 100 % | **115** (32.1 %) |
| Solved ≥ 99 % | **340** (95.0 %) |
| Solved ≥ 96 % | **353** (98.6 %) |
| Average satisfaction | **99.58 %** |
| Largest perfect solve | **354 890 clauses** (Clique Coloring) |
| Fastest >10 K‑clause perfect solve | **22 521 clauses** in **33 ms** (clique\_6\_v40) |

---  

*All timings were measured on the hardware noted in each table, using the default NitroSAT configuration shipped in the repository.* AMD Ryzen 5 5600H with Radeon Graphics (12) @ 4.280GHz single core of this.
