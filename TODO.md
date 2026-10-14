# Charging Effect Simulation

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

___

### QUESTIONS

* Key to Diplomandenraum, part-time employment ending mid November
* Semiconductor Exam on Wednesday instead of Monday?

* Bias, Chuck explanation
* -> "The wafer has no DC path to ground, since an RF-biased chuck is fed through a blocking capacitor and an unbiased wafer floats, so ions and electrons arrive at the open surface with equal currents "

* "Silicon connected to the substrate is the reference at 0 V and stores no charge."
* -> Does that mean charging only occurs on the mask? All deflection comes from charged mask?

* How can I get the energy and angle distribution of the ViennaPS process model?
* Maxwellian flux for electrons? power-cosine?
* Do we track hits of neutral particles (i.e. reflected ions after their initial hit)?
* How does mapping the charge from the previous time step onto the new surface work?
* ad Surface evolution: handled by ViennaPS? what is already present in ViennaPS? What parameters need to be handed over to ViennaPS for surface evolution? flux, energy, impact angle?

### TODO

  * Calculate electron initial velocity using a thermal energy of T_e = 4eV.
  * Launch equal numbers N of ions and electrons. Set macro-weights for both ions and electrons to 1.0 (or maybe leave it out entirely, i think it is not needed). 
  * ad deposited charge: "Each simulated particle represents many real ones, carrying the weight w = ionFlux * A * dt_c/N. ionFlux is the physical ion flux onto the wafer in ions per m2*s (the ionFlux parameter of the ViennaPS model), A is the area of the source plane (in 2D its length times unit depth), dt_c is the time interval that one iteration represents, and N is the number of ions, and of electrons, launched per iteration." 
  * Prevent charge buildup on the silicon. Add a check so that charge hitting silicon drains to 0V and is not stored.
  * Use periodic or reflecting lateral BCs, and a zero-gradient top at least three trench widths above the mask.
  * Replace the hardcoded for loop with a while loop that ends when the change in surface charge between iterations drops below some threshold (i.e. steady-state).
  * After steady state: additional iterations with fixed potential (?) for cycle-averaged ion yields (flux, energy, impact angle)
  * implement ion reflection -> neutralized ion, straight path
  * Include material parameters (permittivity)
  * Solve the Poisson equation in both gas and solid cells, with permittivities in solids. ("Silicon connected to the substrate is the reference at 0 V and stores no charge.")
  * Apply charge of each hit into the dielectric solid cell at the impact point instead of gas cell in front of it







* [ ] **Implement Electron Transport:**

  * [x] Add electron flux

  * [x] Option A: Monte Carlo electrons ($q=-1$, low mass, broad angle, micro-timesteps).

  * [ ] Option B: Analytical geometric view-factor electron model (Eq 5 from paper).

* [ ] **Physical Sheath Velocities:**

  * [x] Replace hardcoded `maxVel` in `Parameters.txt` with sheath-potential-derived initial velocities ($v_0 = \sqrt{2q|E_{ext}|d_{sh} / m}$).

* [ ] **Explicit Charge-Continuity Equation:**

  * [ ] Refactor particle hits into normalized Fluxes ($\Gamma_i, \Gamma_e$).

  * [ ] Update charge using a pseudo-time step ($\Delta t_c$) and Eq. 4 from the literature summary.

* [ ] **Charge Leakage Model:**

  * [ ] Implement surface leakage/decay ($\nu_{leak}$).

  * [ ] Implement bottom-tapering near the conductive substrate to prevent field singularities.

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