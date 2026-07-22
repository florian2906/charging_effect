# Charging Effect Simulation

* Fix leakage issue -> check paraview visualization, error on bottom of trench, neg. potential moving upwards, switching with pos. after each other step??

### DONE

* [x] **Framework Initialization:** Integrated ViennaPS, ViennaLS, and parameter parsing.

* [x] **Geometry Generation:** Built 2D trench geometry, surface mesh (DiskMesh), and volumetric grid (CellSet).

* [x] **Particle Tracing:**

  * [x] Generalized Particle kinematics class (`Particle.hpp`).

  * [x] Embree based ray-tracing integration (`GeometryCollider.hpp`).

  * [x] Power-Cosine angular distribution sampling for directional ion beams.

  * [x] Safe boundary conditions (kill bounds for bottom, sides, and top reflections).

* [x] **Electrostatics:**

  * [x] Mapped 1D surface charge to 3D volumetric ChargeDensity grid.

  * [x] Built Jacobi Poisson Solver for Electric Potential ($\Phi$) (`PoissonSolver.hpp`).

  * [x] Calculated Electric Field gradients ($E_x$, $E_y$).

* [x] **Self-Consistent Loop:**

  * [x] Batched quasi-steady iteration loop (tracing $\rightarrow$ charging $\rightarrow$ field solving).

  * [x] Fixed array bloat bugs in VTK exports.

  * [x] Exporting distinct .vtp and .vtu files per iteration for ParaView.

* [x] **Explicit Charge-Continuity Equation:**

  * [x] Refactor particle hits into normalized Fluxes ($\Gamma_i, \Gamma_e$).

  * [x] Update charge using a pseudo-time step ($\Delta t_c$) and Eq. 4 from the literature summary.

* [x] **Charge Leakage Model:**

  * [x] Implement surface leakage/decay ($\nu_{leak}$).

  * [x] Implement bottom-tapering near the conductive substrate to prevent field singularities.

___

### TODO

* [ ] **Implement Electron Transport:**

  * [x] Add electron flux

  * [x] Option A: Monte Carlo electrons ($q=-1$, low mass, broad angle, micro-timesteps).

  * [ ] Option B: Analytical geometric view-factor electron model (Eq 5 from paper).

* [ ] **Physical Sheath Velocities:**

  * [x] Replace hardcoded `maxVel` in `Parameters.txt` with sheath-potential-derived initial velocities ($v_0 = \sqrt{2q|E_{ext}|d_{sh} / m}$).

* [ ] **Material-Aware Electrostatics:**
  
  * [ ] Upgrade the Poisson solver to recognize different material IDs from the volumetric grid.
  
  * [ ] Apply specific dielectric constants ($\epsilon_r$) for the mask (e.g., 3.9 for SiO2) and gas/vacuum.

* [ ] **Link to reactor scale simulation:**
  
  * [ ] Use ion and electron fluxes from reactor scale simulation instead of randomizing? 
  
  * [ ] Whats the actual ratio of ions/electrons launched?
  
  * [ ] Right now ions are just protons; is that true?

* [ ] **Profile Evolution:**

  * [ ] Feed the steady-state, field-modified fluxes into ViennaPS surface reaction models.

  * [ ] Move level-set boundary to simulate actual physical etching.

* [ ] **Others:**

  * [ ] Fix `WARNING: PointData attempted to access scalar data...` 
  
  * [ ] Consider cleaning up `CharginEffect.cpp`, moving functional blocks into `Helpers.hpp` 
  
  * [ ] Clean up ambigous loop counters
  
  * [ ] Add proper comments, docstrings, update `README.md` 

  * [ ] Benchmark, find bottle neck and improve

  * [ ] GPU parallelization? -> Ray tracing bottle neck?


SF6O2 Plasma Etchcing
surface charge as coverage