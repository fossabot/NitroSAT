# Changelog

All notable changes to NitroSAT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.0.0] - 2026-03-08

### Added

- **NitroSAT v2** - Major rewrite with physics-informed optimizer improvements

### Changed (v1 → v2)

#### 1. Optimizer: NADAM → WAdam (Wasserstein-inform

| Feature | v1 (NADAM) | v2 (WAdam) |
|---------|------------|------------|
| **Algorithm** | Nesterov-accelerated Adam | 1D Wasserstein Flow with Resonance (WFR) |
| **Momentum** | Single `m` vector | Dual: `m_standard` + `m_phase` |
| **Variance** | Single `v` vector | Dual: `v_amp` + `v_phase` |
| **Phase tracking** | None | Tracks gradient direction bias via `prev_phase` |
| **Directional damping** | None | Uses `cos(m_phase)` to cancel oscillating forces |
| **Preconditioning** | Standard Adam | Variance interpolation: `κ * v_amp + (1-κ) * v_phase * v_amp/π²` |

**v2 Code** (`src/c/v2/nitrosatv2.c:232-286`):
```c
// Decompose into Amplitude and Phase-Velocity momentum
opt->m_standard[i] = b1 * opt->m_standard[i] + (1.0 - b1) * grad;
opt->v_amp[i] = b2 * opt->v_amp[i] + (1.0 - b2) * grad * grad;

// Phase momentum tracks directional bias
double phase_cur = (grad < 0) ? PI : 0.0;
double phase_diff = phase_cur - opt->prev_phase[i];
opt->m_phase[i] = b1 * opt->m_phase[i] + (1.0 - b1) * phase_cur;
opt->v_phase[i] = b2 * opt->v_phase[i] + (1.0 - b2) * phase_diff * phase_diff;

// Directional damping - cancels oscillating forces at resonance
double m_wfr = fabs(m_hat) * cos(m_phase_hat);
```

#### 2. Learning Rate Schedule

| Feature | v1 | v2 |
|---------|----|----|
| **Formula** | `A0 * I(t) * ζ-gain` | `A0 * max(0, I(t))` |
| **Negative values** | Capped to `fmax(0, ...)` | Clamped with explicit `if (I < 0.0) I = 0.0` |
| **ζ-gain** | `1 / (1 - progress + 0.02)` | Removed (simpler) |

**v2 Code** (`src/c/v2/nitrosatv2.c:1862-1875`):
```c
static double annealing_lr_v2(double t, double A0) {
    double I = sin(PI/20.0 * t) * exp(-PI/20.0 * t)
             + sin(PI/10.0 * t) * exp(-PI/10.0 * t)
             + sin(PI/9.0  * t) * exp(-PI/9.0  * t)
             + sin((t / (phi * phi))) * exp(-(t / (phi * phi)));
    if (I < 0.0) I = 0.0;
    return A0 * I;
}
```

#### 3. Temperature / β Dynamics

| Feature | v1 | v2 |
|---------|----|----|
| **β schedule** | Fixed `heat_beta` | `compute_beta_eff_v2()` with ζ(s) dynamics |
| **ζ approach** | Direct `1/(s-1)` pole | Self-consistent: `s = 1 + ε * exp(-step/τ) + s_cumulative` |
| **Euler-Maclaurin** | Not used | Uses ζ(s) ≈ `1/(s-1) + γ` approximation |

**v2 Code** (`src/c/v2/nitrosatv2.c:1882-1896`):
```c
static double compute_beta_eff_v2(double beta, int step, int max_steps, ...) {
    double tau = max_steps / 4.0;
    double epsilon = 0.5;
    double s = 1.0 + epsilon * exp(-step / tau) + s_cumulative;
    double zeta_s = 1.0 / (s - 1.0) + GAMMA;
    double zeta_dynamics = (zeta_s + 0.1 * zeta_prime_s) / fmax(fabs(zeta_s), 0.1);
    return beta * (1.0 + zeta_dynamics);
}
```

#### 4. Data Structures

| Feature | v1 | v2 |
|---------|----|----|
| **Unsatisfied clause tracking** | Scan all clauses each time | `unsat_list[]` + `unsat_pos[]` for O(1) updates |
| **Memory** | O(M) full scan per iteration | O(M) build + O(1) incremental updates |
| **CSR validation** | Basic bounds checking | Full validation: clause occurrence count matching |

**v2 Code** (`src/c/v2/nitrosatv2.c:1901-1918`):
```c
static void recompute_sat_counts(NitroSat *ns) {
    ns->unsat_sz = 0;
    for (int c = 0; c < ns->num_clauses; ++c) {
        // ... compute sat_counts[c] ...
        if (sc == 0) {
            ns->unsat_pos[c] = ns->unsat_sz;
            ns->unsat_list[ns->unsat_sz++] = c;
        } else {
            ns->unsat_pos[c] = -1;
        }
    }
}
```

#### 5. Topological Repair

| Feature | v1 | v2 |
|---------|----|----|
| **β₁ calculation** | `β₀ - χ` | `β₀ - χ + unsat_sz` (includes unsatisfied clause count) |
| **Edge counting** | Standard | Adjusted for unsatisfied clause density |

**v2 Code** (`src/c/v2/nitrosatv2.c:672`):
```c
int chi = active_cnt - (int)edge_cnt + ns->unsat_sz;
topo.beta1 = (topo.beta0 - chi > 0) ? (topo.beta0 - chi) : 0;
```

#### 6. Local Search (BAHA Walksat)

| Feature | v1 | v2 |
|---------|----|----|
| **Walksat integration** | Yes, `baha_walksat()` | Removed (pure gradient descent) |
| **Fallback** | Phase 4: walksat for stubborn unsat | Only gradient + topological repair |

#### 7. Performance Features

| Feature | v1 | v2 |
|---------|----|----|
| **Grokking detection** | Yes | Removed |
| **Stagnation monitoring** | Sliding window (50 steps) | Simplified |
| **Nuclear perturbation** | Yes | Removed |

#### 8. Output Format

| Feature | v1 | v2 |
|---------|----|----|
| **Timing breakdown** | `walksat_ms` | Removed (no walksat) |
| **JSON latency fields** | Full phase breakdown | Simplified |

---

## [1.0.0] - 2026-01-15

### Added

- Initial release of NitroSAT (v1)
- O(M) linear time complexity MaxSAT solver
- Continuous relaxation with physics-informed optimization
- Prime-weighted clause learning
- BAHA (Branch-Aware Holonomy Annealing)
- Three-phase finisher: Langevin → Topological Repair → Adelic Saturation
- DRAT proof generation for UNSAT instances
- JSON output with diagnostic information

---

## Comparison Summary

| Aspect | v1 | v2 |
|--------|----|----|
| **Optimizer** | NADAM | WAdam (Wasserstein Flow) |
| **LR Schedule** | ζ-gain modulated | Simple non-negative |
| **Temperature** | Fixed β | ζ(s) self-consistent |
| **Unsat tracking** | Full scan | Incremental list |
| **Walksat** | Yes | No |
| **Complexity** | More features | Streamlined |
| **Use case** | General purpose | Faster on easy cases |

### When to Use Which

- **v1**: For hard instances requiring walksat fallback
- **v2**: For faster execution on easier problems where gradient descent converges

---

*This changelog was generated on 2026-03-08.*
