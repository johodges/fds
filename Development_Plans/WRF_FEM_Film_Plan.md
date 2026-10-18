# Development Plan: WRF Two-Way Coupling, FEM Solid Heat Transfer, and a Surface Liquid-Film Solver

*Baseline: `firemodels/fds` master at commit `9bb47f3` (fork `johodges/fds`). Every `file:line` reference is to that commit.
Effort figures are rough order-of-magnitude estimates for an experienced FDS developer. Refine them after Phase 0.*

This plan covers three features:

- **A.** Two-way coupling with WRF, similar in spirit to WRF-SFIRE.
- **B.** Replacing the solid-phase heat transfer with a finite-element (FEM) or equivalent unstructured representation for complex 3-D conduction.
- **C.** An interfacial thin-film solver for water flowing on surfaces, similar to the FireFOAM surface-film model.

It is based on a survey of the current source (`Source/*.f90`), the four manuals, and the Verification/Validation suites. Sections A.1, B.1 and C.1 record what exists today. The remaining sections propose designs, code touch points, phases, verification and risks.

---

## 1. Executive summary

### A. WRF two-way coupling (rough effort: 14–20 person-months FDS side, plus 3–5 WRF side)

- **Approach.** FDS runs as a microscale LES *nest* inside WRF.
  - WRF → FDS: time-varying 3-D lateral/top boundary data and initial fields.
  - FDS → WRF: fire heat, water-vapour and smoke *sources*, binned onto the WRF grid. They are injected through WRF-Fire's tendency path, which is the SFIRE analogue.
  - Staging: offline one-way first, then online two-way (MPMD).
- **Main reuse.**
  - OPEN-boundary machinery, DFSEM, and the `PATCH_VELOCITY_FLUX` forcing pattern.
  - CSVF initial fields, `UPDATE_HRR` integrals, and the external-control heartbeat.
  - The unused `ORIGIN_LAT/ORIGIN_LON/NORTH_BEARING` inputs.
- **Main obstacles.**
  - Wind is a single U(z) profile per mesh. Inflow T and species are uniform.
  - `MPI_COMM_WORLD` is hard-coded 307 times.
  - There is no geolocation.
  - DT is not clipped to exchange times.

### B. FEM solid heat transfer (rough effort: 18–28 person-months)

- **Approach.** A new unstructured 3-D conduction and pyrolysis solver that **replaces HT3D**.
  - P1 tetrahedra with a lumped mass matrix, which is equivalent to vertex-centred CVFEM.
  - Implicit time stepping, with HYPRE PCG/AMG.
  - The 1-D solver stays as the default for 1-D and thin surfaces.
- **Main reuse.**
  - Tetrahedra storage (`GEOMETRY%VOLUS/MATLS`), `BOX2TETRA`, and `.bingeom`.
  - `PROC_MATL` tables, the `PYROLYSIS` kinetics, and `HYPRE_INTERFACE`.
  - The HT3D verification suite, which includes an ANSYS FEM I-beam benchmark.
- **Main obstacles.**
  - HT3D is a direction-split set of 1-D lines with replicated state.
  - The data layout is face-centric, and the code assumes Cartesian `IOR` orientations.
  - `PYROLYSIS` reads gas arrays directly.
  - The solid partition is tied to the gas meshes.

### C. Surface liquid film (rough effort: 14–22 person-months)

- **Approach.** FireFOAM-style depth-averaged film mass, momentum and energy equations.
  - Finite volumes on a dedicated **surface mesh**: wall-cell quads for OBSTs, GEOM triangles for complex geometry.
  - The film is an extra *implicit node* in the solid's surface energy balance.
  - Coupled to Lagrangian sprays (absorb, drip, separate, splash) and to the gas through the existing wall mass-flux path.
- **Main reuse.**
  - Surface-droplet hooks in `MOVE_PARTICLES` and `CFACE_SEARCH`.
  - Droplet evaporation correlations and the wall-function mass transfer in `CONDENSATION_EVAPORATION`.
  - The `CALCULATE_ZZ_F` mass-flux path and `EDGE_FACES` triangle adjacency.
  - The back-wall MPI exchange pattern.
- **Main obstacles.**
  - Wall cells have no surface connectivity, and CFACEs have no adjacency.
  - A single surface can mix WALL cells and CFACEs.
  - Today's droplet–wall coupling is explicit and lagged.

### 0. Shared foundations (rough effort: 3–5 person-months)

A communicator abstraction, a surface-mesh layer, a sparse-solver wrapper, and an explicit solid–gas exchange contract.

### Recommended order

1. **Phase 0 first.** It is small and unblocks all three tracks.
2. **Then run the three tracks in parallel.** A1 (offline WRF → FDS) and B1 (serial FEM conduction) deliver value earliest.
3. **Film work.** C starts on Cartesian walls with 1-D solids and moves onto FEM solid surfaces once B5 lands.

---

## 2. Phase 0: Shared foundations

### 0.1 An FDS communicator

**Current state.**
- `MPI_COMM_WORLD` appears 307 times in 11 files: main 152, ccib 74, pres 26, geom 26, dump 7, fire 6, hvac 5, read 4, func 3, turb 2, radi 2.
- FDS owns MPI start-up and shut-down:
  - `MPI_INIT_THREAD` at `main.f90:89`.
  - `MPI_FINALIZE` followed by `STOP` in `END_FDS` (`main.f90:1794-1869`).
  - `MPI_ABORT` on WORLD at `main.f90:3987` and `ccib.f90:5124`.
- Communicators derived from WORLD include `MPI_COMM_NEIGHBORS` (`main.f90:1238-1260`) and `ZSL_COMM` (`pres.f90:4470`). The MKL cluster solver is also handed WORLD (`pres.f90:3361` ff.).

**Change.**
- Add `TYPE(MPI_COMM) :: MPI_COMM_FDS` to `cons.f90`, next to `MPI_COMM_NEIGHBORS` (`cons.f90:389`).
  - Standalone: it equals `MPI_COMM_WORLD`.
  - Coupled: it is FDS's colour of an `MPI_COMM_SPLIT`.
- Take `MY_RANK` and `N_MPI_PROCESSES` from it, and replace every use of WORLD. The replacement is mechanical because the code uses `MPI_F08`.
- Skip `MPI_INIT` if `MPI_INITIALIZED` is already true. In coupled mode, let the coupler finalize.

**Acceptance.** The whole verification suite is bit-identical.

**Needed by.** A2 (MPMD with WRF), B3 and C3 (sub-communicators for the solid and film ranks).

### 0.2 Surface-mesh layer (new module, e.g. `smsh.f90`)

A conformal polygonal surface mesh with edge adjacency. It is independent of the gas grid and carries conservative maps to the gas-side boundary faces.

**Contents.**
- Facets: area, unit normal, centroid, vertices, owner mesh/rank.
- Edges: two facets, fold angle, edge vector.
- Maps from facet to `WALL(IW)` or `CFACE(ICF)`, with area weights.

**Builders.**

1. **Cartesian OBST surfaces.** One quad per solid WALL cell. Neighbours across the four edges are derived from `CELL%WALL_INDEX(-3:3)` (`type.f90:2174-2186`), for a gas cell G, face normal n and tangential direction d:
   - Planar neighbour: `CELL(G+d)%WALL_INDEX(-n)`.
   - Concave corner: `CELL(G)%WALL_INDEX(d)`.
   - Convex edge: `CELL(G+d-n)%WALL_INDEX(-d)`.

   This is the topology `MOVE_PARTICLES` follows implicitly today (`part.f90:2515-2553`). Cross-mesh neighbours use `FIND_WALL_INDEX` (`func.f90:5501`) and the remote `CELL` data from `EXCHANGE_GEOMETRY_INFO` (`func.f90:5377`).
2. **GEOM surfaces.** Use `GEOMETRY%FACES` with `FACE_EDGES`/`EDGE_FACES` (`type.f90:1191`; built by `GET_GEOM_EDGES`, `geom.f90:946`).
   - Closed manifold geometry is already enforced (`geom.f90:876-893`), so every edge has two triangles.
   - Add optional refinement down to about one gas cell.
   - Map facets to CFACEs through `CUT_FACE%BODTRI` (`type.f90:1305`) with area weights. The map is many-to-many because cut cells can be merged or blocked.
3. **FEM boundary triangles** (Part B).

**Halo exchange.** A new `MESH_EXCHANGE` code with persistent requests. Model it on codes 66/67 (`main.f90:3341-3372`, `3726-3760`) and `INITIALIZE_BACK_WALL_EXCHANGE` (`main.f90:2204-2460`).

**v1 is static.** Stop with an error if a filmed or FEM surface belongs to a removable OBST. Dynamic rebuild later hooks into `REASSIGN_WALL_CELLS` (`init.f90:4909`).

**Needed by.** C (the film domain) and B5 (the FEM–CFACE interface).

### 0.3 Sparse linear-solver wrapper

- Generalize `HYPRE_INTERFACE` (`imkl.f90:200-543`) into a small API:
  - IJ assembly, including off-process entries.
  - PCG (SPD) and GMRES, with BoomerAMG and setup reuse.
  - A communicator argument.
  - **Runtime** tolerances and iteration limits. Today they are compile-time `PARAMETER`s at `imkl.f90:204-229`.
- Keep PARDISO (`imkl.f90:7-96`) as an optional single-rank direct solver.
- Copy the global unknown numbering from `GET_MATRIX_INDEXES_H` (`pres.f90:5406-5513`).
- HYPRE is already fetched by default in `CMakeLists.txt` (`USE_HYPRE` ON), so this adds no new dependency.

### 0.4 Solid–gas exchange contract

Make explicit the per-face contract that the 1-D solver already satisfies. The FEM and film solvers can then plug in without touching gas-phase code.

**Inputs from the gas.** `B1%TMP_G`, `HEAT_TRANS_COEF`, `Q_RAD_IN`, `ZZ_G`, `Q_CONDENSE`, `B2%LP_CPUA`.

**Outputs to the gas.**
- `B1%TMP_F`.
- `Q_CON_F`, time-averaged over sub-steps (`wall.f90:2920-2928`, `2968-2969`).
- `Q_RAD_OUT`, `EMISSIVITY`, `M_DOT_G_PP_ADJUST/ACTUAL`, `Q_DOT_G_PP`, `M_DOT_PART_ACTUAL`, `T_IGN`, `BURNAWAY`.

They are consumed by `CALCULATE_ZZ_F` (`wall.f90:1036`), `divg.f90:532-554`, `ccib.f90:10948-11119` and `radi.f90:4220-4250`.

**Refactor of the pyrolysis kinetics.** Split the pointwise kinetics in `PYROLYSIS` (`wall.f90:3197-3607`) from its direct gas-array reads (`ZZ`, `TMP`, `PBAR`, `RSUM`, `RDX` at `IIG/JJG/KKG`). The result is a pure routine of (T, ρ_i, local gas state) that both the 1-D and FEM solvers call.

### 0.5 Ground rules for all three tracks

- **Opt-in only.** Every feature is off by default, and default inputs stay bit-identical (Firebot regression).
- **Each phase ships with:**
  - Verification cases registered in `Verification/FDS_Cases.sh`.
  - User Guide and Technical Reference Guide sections.
  - Validation cases where data exist.
- **New files.**
  - Follow the 4-letter naming used by the existing sources: `cplr.f90` (coupler), `sfem.f90` (solid FEM), `film.f90`, `smsh.f90`.
  - Add each one to both the explicit list in `CMakeLists.txt` and `obj_mpi` in `Build/makefile`, in dependency order.

---

## 3. Part A: Two-way coupling with WRF

### A.1 What exists today

**Wind and atmosphere.**
- `&WIND` is read at `read.f90:2091-2331`. It builds **one vertical profile per mesh**:
  - The profile is stored in `U_WIND/V_WIND/W_WIND(0:KBP1)` (`mesh.f90:100-102`).
  - It is filled by `COMPUTE_WIND_COMPONENTS` (`func.f90:2827-2840`) from `GET_WIND_AT_HEIGHT` (`func.f90:2679-2705`).
- Physics already available:
  - Monin–Obukhov profiles (`func.f90:2539-2574`), pressure-gradient force and `FORCE_VECTOR`, geostrophic wind, and Coriolis (`velo.f90:860-1031`).

**OPEN boundaries.**
- Tangential ghost velocity = `U_WIND(K)` plus SEM eddies (`velo.f90:2112-2198`).
- Normal velocity is set only indirectly, through the Dirichlet H condition (`pres.f90:147-224`).
- Inflow temperature is `TMP_0(z)` or a single `TMP_EXTERIOR` per VENT. Inflow species are the uniform ambient values (`wall.f90:668-721`).

**Background state.**
- `TMP_0/P_0/RHO_0(z)` are static. They come from a hydrostatic integration at start-up (`init.f90:438-487`).
- `HUMIDITY` is a single value folded into the `AIR` composition (`read.f90:3169`, `3828-3830`).

**Other building blocks.**
- **Synthetic turbulence.** DFSEM (`turb.f90:1623-1909`). Eddies are advected by one height-mean wind (`turb.f90:1726-1738`).
- **In-domain relaxation template.** `PATCH_VELOCITY_FLUX` (`velo.f90:1065-1200`) drives the velocity toward a target with IBM-like forcing.
- **Initial 3-D fields.** The CSVF route (`init.f90:5322-5539`).
- **External control.**
  - Rank 0 reads scalar RAMP/CTRL values from a CSV every `DT_EXTERNAL` (`func.f90:545-622`), with an optional heartbeat.
  - A ramp marked external loses its z-dependence (`func.f90:826-828`).

**Geolocation.**
- `ORIGIN_LAT`, `ORIGIN_LON` and `NORTH_BEARING` (`cons.f90:361-363`) are only echoed to the `.out` file (`dump.f90:2817-2821`).
- There is no projection. The code assumes x = east and y = north.

**Outputs usable for feedback.**
- Per-cell `Q`, `QR` and `M_DOT_PPP`.
- Per-wall `Q_CON_F` and `M_DOT_G_PP_ADJUST`.
- The `UPDATE_HRR` integrals (`dump.f90:10455-10641`), which are a template for per-column aggregation.

**History.**
- A one-way WRF mean-forcing experiment existed in 2014 (`MEAN_FORCING`, `RAMP_U0`, "z profile ramps for mean forcing, WRF coupling"; `Utilities/Misc/GoogleCode_Commit_Log.txt:2560-2568`). It was later removed.
- `Verification/Atmospheric_Effects/wind_example_5/10/32.fds` still carry "one-way WRF coupling" titles.
- **Ask the NIST developers why it was dropped before re-introducing nudging.** The commit log hints at MPI problems with nudging.

### A.2 Target architecture

```
 ┌──────────── WRF parent (η coordinates, map projection, possibly LES-resolving) ─────────────┐
 │ solve_em advances [t_n, t_n+1] using the FDS fire sources from the previous interval        │
 └──────────────┬──────────────────────────────────────────────────────────────▲──────────────┘
                │ u, v, w, θ, q_v, p, PBL diagnostics                         │ heat, H2O, smoke sources
                │ on FDS boundary + relaxation points, at t_n+1               │ on WRF cells, averaged over [t_n, t_n+1]
 ┌──────────────▼──── coupler layer cplr.f90 (backends: OFFLINE | OASIS3-MCT | MPI) ───────────┴──────────────┐
 │ interpolation/rotation weights precomputed by Utilities/Python/wrf2fds (projection, AGL, terrain blend)      │
 └──────────────┬──────────────────────────────────────────────────────────────▲──────────────────────────────┘
                ▼                                                              │
 ┌───────────────────────── FDS nest (Cartesian low-Mach LES) ──────────────────┴──────────────────────────────┐
 │ per-face boundary arrays at t_n and t_n+1, interpolated in time → OPEN BCs for velocity, T, Y_H2O           │
 │ inflow turbulence (DFSEM or cell perturbation) + optional lateral velocity relaxation zone                   │
 │ DT clipped to land on t_n+1; source accumulators (q''' + QR, wall Q_CON_F, species ṁ''')                     │
 └──────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### A.3 Key design decisions

**1. FDS is a nest, not an embedded fire module.**
- WRF-SFIRE embeds a semi-empirical level-set spread model on a refined fire grid inside WRF. Because WRF has no flux boundary conditions, the fire's sensible and latent heat is inserted as forcing terms in a layer above the surface, with an assumed exponential decay in height (Mandel et al. 2011).
- FDS's value is its resolved fire and plume physics. So FDS should receive WRF boundary data and return sources.
- An **SFIRE-equivalent degenerate mode** is still worth building as a stepping stone (A2a). In it, `LEVEL_SET_MODE=1` is driven by WRF winds and returns surface fluxes. It isolates the coupling infrastructure from LES issues and can be compared directly with WRF-SFIRE.

**2. What WRF sends.**
- Sent on FDS lateral and top boundary points and on relaxation-zone cells:
  - u, v, w (earth-relative, then rotated into FDS axes).
  - T, from θ (θ_m in WRF v4) and p = P + PB.
  - Y_H2O = q_v/(1+q_v).
  - PBL diagnostics (u\*, TKE or PBL height, surface heat flux) for the inflow turbulence generator.
  - Surface pressure.
- Full 3-D initial fields are sent once.

**3. What FDS returns (feedback mode `SOURCES`, the default).** For each WRF cell that overlaps the FDS domain, interval-averaged:
- **Sensible heating rate**:
  - the sum of (q''' + QR)·dV over the FDS gas cells inside it;
  - plus wall convective heat to the gas (−Q_CON_F·A in FDS's sign convention);
  - minus droplet-evaporation sinks.

  WRF turns this into a θ tendency.
- **Water-vapour source** (combustion, fuel moisture, suppression water), as a q_v tendency.
- **Smoke/soot**, and optionally CO and CO₂ for WRF-Chem, as tracer tendencies.

This generalizes SFIRE's assumed exponential profile with FDS's resolved 3-D distribution. WRF still carries the plume beyond the FDS box. The fire-induced far-field indraft then returns to FDS through its lateral boundaries, which closes the two-way loop.

A **2-D variant** sends column integrals only. These map directly onto WRF-Fire's ground and canopy heat/moisture flux arrays, so a first version can reuse WRF-Fire's existing injection code unchanged.

**Optional feedback mode `STATE` (later).** WRF prognostic fields in the overlap are relaxed toward box-filtered FDS fields, as in WRF's own two-way nest feedback. It is mutually exclusive with `SOURCES`, which avoids double counting.

**4. Coordinates and geolocation.**
- Give `ORIGIN_LAT`, `ORIGIN_LON` and `NORTH_BEARING` real meaning, and add `ORIGIN_ELEVATION` for the vertical datum.
- The preprocessor maps WRF projection coordinates (Lambert, Mercator, polar stereographic, lat-lon) to local x, y.
- Winds are rotated grid-relative → earth-relative (WRF `COSALPHA/SINALPHA`) → FDS axes (`NORTH_BEARING`).
- All interpolation and rotation weights are precomputed and stored as a SCRIP-style weights file. At runtime FDS only does sparse dot-products.

**5. Vertical interpolation and terrain mismatch.**
- WRF heights come from (PH+PHB)/g on w-levels. WRF terrain (`HGT`) is much smoother than the FDS DEM terrain.
- Interpolate in height **above ground** at each boundary column.
- Offer a terrain-blending strip along the FDS lateral boundaries, where FDS terrain relaxes to WRF terrain. Do this in preprocessing, next to qgis2fds.
- For FDS points below WRF's lowest level, extrapolate with a Monin–Obukhov profile built from WRF's u\* and L.

**6. Thermodynamics.**
- Set `P_INF` and `TMP_0(z)` from the WRF footprint-mean column at start-up, using the existing `RAMP_TMP0_Z` route.
- Hold the background frozen in v1. This is a documented limitation for multi-hour runs. A time-varying P_0 needs the dP̄/dt term in the divergence and comes later.
- Variable humidity needs an **explicitly tracked WATER VAPOR species over a dry-AIR background**, because `HUMIDITY` is folded into AIR today.

**7. FDS boundary conditions.**
- Replace the 1-D `U_WIND(K)` with per-external-wall-cell target arrays at two time levels: `U_BC`, `V_BC`, `W_BC`, `TMP_BC`, `ZZ_BC(:)`. Interpolate linearly in time.
- Keep the existing OPEN formulation unchanged except that it now works per face:
  - Dirichlet H for the normal velocity.
  - Ghost-cell tangential velocity.
  - Upwind selection of T and Y.
- If the soft OPEN coupling drifts in strong cross-flow, evaluate a specified-normal-velocity inflow with a global mass-flux correction on the inflow faces, as PALM's offline nesting does.

**8. Inflow turbulence.**
- Mesoscale WRF fields carry no resolved turbulence. Two options:
  - DFSEM with per-face statistics scaled from WRF PBL diagnostics. Today it uses one mean wind and one intensity (`turb.f90:1726-1738`, `2021-2023`).
  - The stochastic cell-perturbation method of Muñoz-Esparza et al., applied to θ in a boundary strip. It is cheap and widely used for WRF → LES nesting.
- When the parent is WRF-LES, fluctuations are resolved, and the exchange interval should be ≤ the WRF time step.

**9. Lateral relaxation zone (optional).**
- Davies-type relaxation of velocity toward WRF over N cells inside the boundary.
- Implement it like `PATCH_VELOCITY_FLUX` and call it next to `DIRECT_FORCE` (`velo.f90:839-841`).
- It absorbs plume-induced reflections. It covers velocity only in v1, because scalar relaxation would have to pass through the divergence and the species sources.

**10. Time synchronization: sequential, staggered, explicit.**
```
for each coupling interval [t_n, t_n+1]   (Δc = integer × WRF dt):
    WRF advances t_n → t_n+1 using FDS sources averaged over [t_n-1, t_n]
    WRF → FDS: boundary state at t_n+1
    FDS advances t_n → t_n+1, BCs interpolated in time, DT clipped to land exactly on t_n+1
    FDS → WRF: sources accumulated over [t_n, t_n+1]
```
- The feedback lags by one interval, much like SFIRE's explicit coupling.
- FDS dominates the cost, so give WRF few ranks.
- A lag-2 variant can overlap the two codes' execution later.
- Before ignition, spin up FDS turbulence for 10–30 minutes with feedback off.

**11. Coupler backends.**
- An abstract `cplr.f90` API: `CPL_INIT`, `CPL_RECV_BOUNDARY`, `CPL_SEND_SOURCES`, `CPL_FINALIZE`.
- Backends:
  - **(a) OFFLINE.** Pre-generated boundary files; one-way.
  - **(b) OASIS3-MCT.** Parallel regridding with precomputed weights and MPMD launch.
    - WRF already contains an OASIS3-MCT interface for ocean coupling (`frame/module_cpl*.F` in recent releases). **Confirm this for the target WRF version.**
  - **(c) Plain MPI intercommunicator.** The fallback.
- preCICE and ESMF/NUOPC are alternatives, but they are heavier for 3-D volume data.

### A.4 FDS-side work by file

| File | Change |
|---|---|
| `cons.f90`, all MPI users | `MPI_COMM_FDS` (Phase 0.1). |
| `read.f90` | New `&CPLR` namelist: `MODE` (`OFFLINE`/`ONLINE`), `BOUNDARY_FILE`, `WEIGHTS_FILE`, `DT_COUPLE`, `FEEDBACK` (`NONE`/`SOURCES`/`STATE`), `RELAXATION_ZONE_CELLS`, `INFLOW_TURBULENCE` (`SEM`/`CELL_PERTURBATION`/`NONE`), `ORIGIN_ELEVATION`. Validate the existing MISC `ORIGIN_*`/`NORTH_BEARING`. Require an explicit WATER VAPOR species when coupled. |
| `mesh.f90`, `init.f90` | Per-external-wall-cell boundary arrays at two time levels. Allocate them like `U_GHOST` (`init.f90:972-986`), after `INITIALIZE_ATMOSPHERE` and `INITIALIZE_WALL_ARRAY` (`main.f90:182-185`). Initial fields go through the CSVF path (`init.f90:5322-5539`) or a new binary reader. |
| `func.f90` | Per-face path through `COMPUTE_WIND_COMPONENTS`/`GET_WIND_AT_HEIGHT` (`func.f90:2679-2840`). Time interpolation helper. |
| `velo.f90` | `VELOCITY_BC` OPEN_AND_WIND block (`velo.f90:2112-2198`) uses per-face targets. Relaxation-zone forcing next to `velo.f90:839-841`. DT clip to the coupling time in `CHECK_STABILITY` (`velo.f90:3028-3208`). |
| `pres.f90` | Per-face H0 target velocity (`pres.f90:172-183`). Cut-cell equivalent (`ccib.f90:2540-2542`). |
| `wall.f90` | Inflow T and ZZ in `SURFACE_HEAT_TRANSFER` INFLOW_OUTFLOW (`wall.f90:682-706`). |
| `ccib.f90` | `GET_OPENBC_TANGENTIAL_CUTFACE_VEL` (`ccib.f90:410-444`). |
| `turb.f90` | Per-face SEM advection velocity and intensity (`turb.f90:1726-1738`, `2021-2023`). New cell-perturbation routine. |
| `main.f90` | Coupler init after `READ_DATA` and `MPI_INITIALIZATION_CHORES(3)` (`main.f90:119`, `168`). Receive next to the external-file block (`main.f90:654-661`), outside `CHANGE_TIME_STEP_LOOP` so a DT cut does not repeat it. Clip DT next to the T_END clip (`main.f90:620`). Send after `EXCHANGE_GLOBAL_OUTPUTS`/`UPDATE_CONTROLS` (`main.f90:1037-1041`). Finalize in `END_FDS`. |
| `dump.f90` | Source aggregation modelled on `UPDATE_HRR` (`dump.f90:10455-10641`). A new `CHID_cpl.csv` energy and mass budget of the exchanged fluxes. Coupler state in `DUMP_RESTART`/`READ_RESTART` (`dump.f90:3575`, `3762`). |
| `vege.f90` | SFIRE-equivalent mode: spatially varying WRF wind for `LEVEL_SET_MODE=1` (`vege.f90:356-365`), and surface heat/moisture flux from fuel consumption for the return path. |
| `cplr.f90` (new) | Backend-neutral API, sparse weight application, time interpolation, accumulators. |
| `Utilities/Python/wrf2fds` (new) | Read wrfout (netCDF), compute projection/rotation/AGL weights and the terrain blend, write FDS boundary files, initial fields and the weights file. It keeps FDS itself free of netCDF. |

### A.5 WRF-side work (in a WRF fork, outside this repository)

- **Registry.** Namelist options (e.g. `fds_coupling`, `fds_coupling_interval`, feedback mode) and state arrays for the received sources.
- **Communicator.** Run WRF on its colour of the split. WRF's DM layer accepts an externally supplied communicator, which is how MPMD-coupled systems such as COAWST run it; confirm the mechanism for the target version.
- **Send.** At each coupling time, sample the WRF fields on FDS boundary points; the coupler gathers them from the WRF tiles.
- **Receive.** Add FDS sources as θ, q_v and tracer tendencies.
  - Reuse WRF-Fire's fire-tendency arrays and injection, with WRF-Fire's own spread model switched off.
  - For the 2-D variant, fill WRF-Fire's ground/canopy flux arrays and keep its `fire_ext_grnd`/`fire_ext_crwn` extinction-depth distribution.
  - The array names differ between NCAR WRF (WRF-Fire in `phys/module_fr_fire_*.F`) and OpenWFM WRF-SFIRE. **Pick the target first.**
- **Restart.** Align WRF and FDS restart times.

### A.6 Phases

**A1. Offline one-way, WRF → FDS (5–7 person-months).**
- Build the `wrf2fds` preprocessor and the FDS per-face boundary arrays, reader, time interpolation and DT clip.
- Add geolocation, explicit water vapour, per-face SEM or cell perturbation, and the relaxation zone.
- *Exit:* uniform-profile boundary files reproduce the existing `MO_velocity_profile_*` and `atmospheric_boundary_layer_*` results. Mapping reproduces linear fields exactly. The mass budget closes.

**A2a. Online infrastructure and SFIRE-equivalent mode (3–4 person-months FDS, 2–3 WRF).**
- `MPI_COMM_FDS`, the OASIS3-MCT or MPI backend, and MPMD launch scripts.
- FDS `LEVEL_SET_MODE=1` driven by WRF winds and returning surface fluxes.
- *Exit:* on WRF's idealized fire test, compare against WRF-SFIRE: spread rate, plume and heat-flux totals.

**A2b. Full LES nest with `SOURCES` feedback (3–5 person-months FDS, 1–2 WRF).**
- Energy-budget closure: the heat returned must equal FDS HRR minus radiation lost through the boundaries minus conduction into solids, read from the `.hrr` columns.
- Restart of both codes.

**A3. Validation and hardening (3–4 person-months).**
- FireFlux I/II grass fires (used for WRF-SFIRE evaluation), and RxCADRE if the data can be obtained.
- Sensitivity to coupling interval, FDS domain size and feedback mode.
- The `STATE` feedback option.
- Optional: time-varying background pressure, and dead-fuel moisture driven by WRF RH. Fuel moisture is currently static and applied at start-up only (`vege.f90:56`).

### A.7 Verification and validation

**Verification.**
- Weights reproduce constant and linear fields exactly.
- Binned sources conserve totals between FDS and WRF.
- Restart reproduces a continuous run.
- Results are invariant to the FDS MPI process count.

**Regression.** Feeding the existing `&WIND` profiles through the new file path reproduces:
- `Verification/Atmospheric_Effects/*`
- `Verification/Turbulence/sem_open_wind*`

**Validation, flow only.**
- `Validation/Askervein_Hill` driven by WRF instead of M-O inflow.
- `Validation/Schoenberg_Ekman_Layers` against WRF-LES.

**Validation, coupled fire.** FireFlux; CSIRO grassland cases (`Validation/CSIRO_Grassland_Fires`) with WRF-supplied winds.

### A.8 Risks and mitigations

- **Soft OPEN coupling drifts in strong or rotating cross-flow.**
  - Mitigate with the relaxation zone, or with specified-velocity inflow plus mass correction.
- **WRF and FDS terrain disagree at the lateral boundaries.**
  - Mitigate with AGL interpolation and the blending strip.
- **Inflow turbulence needs a fetch to develop, which inflates domain size.**
  - Mitigate with cell perturbation and by running WRF at LES resolution.
- **Double counting.**
  - Keep `SOURCES` and `STATE` mutually exclusive, and check the budget every interval.
- **Low-Mach vs. compressible formulations.** WRF pressure perturbations enter FDS only through the velocities.
  - Document this, and add a WRF-derived mean pressure-gradient force through `RAMP_PGF_T` if needed.
- **WRF fork maintenance.**
  - Keep the WRF-side changes small: one hook to send, one to receive, reusing the WRF-Fire paths.
- **Scale.** km-scale domains at ~1 m are 10⁹⁺ cells.
  - `SOURCES` feedback works for any FDS domain size relative to the WRF grid.

---

## 4. Part B: FEM (or equivalent) solid-phase heat transfer

### B.1 What exists today

**The 1-D solver.**
- `SOLID_HEAT_TRANSFER` (`wall.f90:1810-3156`) solves one 1-D problem per face along the face normal. Faces are WALL, CFACE, THIN_WALL or particle.
- Numerics:
  - Crank–Nicolson with a Thomas solve (`wall.f90:2854-2910`).
  - A ghost-node surface BC with linearized εσT⁴, and lagged h and q_rad,in (`wall.f90:2867-2894`).
  - Adaptive sub-stepping (`wall.f90:2307-2336`).
- Geometry and properties:
  - Cartesian, cylindrical or spherical geometry.
  - Multi-layer, multi-component.
  - k(T) and c(T) from 1 K tables built in `PROC_MATL` (`read.f90:7524-7869`).
  - Lagrangian shrink/swell with conservative remap (`wall.f90:2366-2815`).
- It runs every `WALL_INCREMENT` gas steps (default 2), in the corrector (`wall.f90:89-99`).

**HT3D is not a 3-D solve.**
- Each exposed face of an HT3D OBST owns a 1-D line through the full depth of the solid.
- Each solver call advances only the lines in one axis direction (`wall.f90:1875`, `1910`), over 3·DT_BC.
- `HT3D_TEMPERATURE_EXCHANGE` (`wall.f90:3615-3725`) then adds the temperature increments to the overlapping nodes of the other two directions. The weights are static and volume-overlap based, built in `INITIALIZE_HT3D_WALL_CELLS` (`init.f90:1831-2284`).
- Consequences:
  - State is replicated per direction and per side.
  - Densities and compositions are never reconciled between lines.
  - Remeshing invalidates the overlap map.
  - Conductivity is isotropic only.

**Documented HT3D limitations** (`FDS_User_Guide.tex:2537-2549`):
- No exterior boundaries.
- No lateral conduction where a 1-D solid contacts a 3-D one.
- No `BURN_AWAY`, and no significant shrink or swell.
- No create/remove.

**GEOM solids** get 1-D through-thickness conduction only.

**What can be reused.**
- `GEOMETRY_TYPE` already stores tetrahedra (`VOLUS`) and `MATLS` (`type.f90:1181-1190`). The `.bingeom` format carries them (`geom.f90:1317-1343`).
- Boundary triangles are extracted from tetrahedra (`geom.f90:2334` ff.). The code carries a comment anticipating this use: "When Tetrahedrons are actually used for heat transfer and pyrolysis this will be needed" (`geom.f90:2322`).
- `BOX2TETRA` (`geom.f90:3266`).

**Gaps in GEOM input.** `VOLUS` and `MATL_ID` cannot be given in `&GEOM` (`geom.f90:1147-1152`), and `MATL_ID` is reset at `geom.f90:3217`.

**Linear algebra.** HYPRE IJ with PCG/BoomerAMG (`imkl.f90:200-543`), and the global numbering in `pres.f90` `GLOBMAT_SOLVER`.

### B.2 Recommendation: replace HT3D, keep the 1-D solver

**Recommendation.** Replace HT3D (and the `VARIABLE_THICKNESS` path) with an unstructured 3-D solver. Keep the 1-D solver as the default for surfaces where 1-D conduction is adequate.

**Why not replace everything.**
- The 1-D solver is fast and embarrassingly parallel.
- It has validated physics paths with no 3-D counterpart: liquid-fuel B-number evaporation, SPyro, TGA, shrink/swell remeshing, thermally thick particles.
- It underpins most of the Validation suite.
- A 1-D column is just a degenerate FEM mesh, so replacing it gains nothing and risks every pyrolysis validation case.

**Migration.**
- Existing HT3D inputs (`SURF HT3D=T` plus `OBST MATL_ID`) should run unchanged on the new solver, so HT3D users move over transparently.
- Keep the legacy HT3D behind `MISC HT3D_SOLVER='LEGACY'` until the FEM solver matches it on the full HT3D verification suite. Then deprecate it.

### B.3 Numerical formulation

**Equation.** ρc(T) ∂T/∂t = ∇·(**k**(T)∇T) + q̇''' in each solid region.
- q̇''' is the sum of reaction heat, the user internal source (`HEAT_SOURCE` with its ramp) and, later, absorbed in-depth radiation.
- **k** is a scalar or a symmetric tensor per element. Anisotropy (wood grain, laminates) is a new capability; HT3D has none.

**Space discretization.** Continuous P1 Galerkin on tetrahedra, with a lumped (diagonal) mass matrix.
- For P1 simplices, the Galerkin stiffness matrix equals the vertex-centred finite-volume (median-dual, CVFEM) diffusion operator.
- So this choice *is* the "equivalent" conservative FV method. Nodal unknowns carry dual-cell volumes, which gives local mass and energy conservation for pyrolysis bookkeeping.
- **Mesh quality.** With lumped mass, the discrete maximum principle needs non-obtuse tetrahedra.
  - Orthoscheme (Kuhn) subdivision of the (graded) box grid used for OBSTs gives non-obtuse tetrahedra.
  - For imported GEOM meshes, mesh-quality checks warn on obtuse dihedral angles.

**Time integration.** A θ-scheme: (M_L/Δt + θA)Tⁿ⁺¹ = (M_L/Δt − (1−θ)A)Tⁿ + F.
- **Default θ = 1** (backward Euler; L-stable, no oscillation at large Fourier numbers).
- Crank–Nicolson as an option, for one-to-one comparison with the 1-D solver.
- BDF2 later.

**Boundary condition.** A Robin condition per boundary facet: −k∇T·n = h(T − T_g) + εσT⁴ − q_rad,in − q_liquid.
- εσT⁴ is linearized as εσ(4T*³T − 3T*⁴), the same trick as the 1-D solver.
- h, T_g and q_rad,in are lagged, as in 1-D.
- Dirichlet (`TMP_FRONT`) and insulated options are supported.
- Exterior-boundary facets use the SURF backing condition, which removes HT3D limitation 1.

**Nonlinearity.**
- Picard iteration (1–3 sweeps, with a convergence test) on k(T) and ρc(T), evaluated from the `PROC_MATL` tables. Interpolate between table entries; the 1-D solver truncates to `INT(T)`.
- An enthalpy-consistent update keeps energy conservation with T-dependent c_p and reactions. The `ht3d_energy_conservation*` cases are the gate.

**Linear solve.**
- The system is SPD because the Robin and mass terms are symmetric and the coefficients are lagged. Solve it with PCG + BoomerAMG through the Phase 0.3 wrapper.
- Reuse the AMG setup across steps while coefficients change slowly.

### B.4 Meshes and input

**1. OBST solids (`SURF HT3D=T`).**
- Generate a box sub-grid that **conforms to the gas wall-cell faces** on exposed surfaces: each wall face equals 2 boundary triangles, or 4 with refinement.
- Grade the sub-grid inward using the existing `CELL_SIZE`, `CELL_SIZE_FACTOR` and `STRETCH_FACTOR` semantics.
- Split each box into tetrahedra.
- Connected OBSTs, including different `MATL_ID`s, are meshed conformally as one region. Material interfaces then have perfect contact.

**2. GEOM solids.**
- Tetrahedra come from a `.bingeom` file (`VOLUS`) or are generated by a new preprocessing utility, `Utilities/Python/fds_tetmesh`, which wraps gmsh or TetGen. It reads STL/OBJ/STEP and writes `.bingeom` with tetrahedra and material tags.
- Enable `MATL_ID` and per-tetrahedron materials in `&GEOM` (remove the reset at `geom.f90:3217`). Restore the orientation and volume constraint flagged at `geom.f90:2322`.
- Do **not** embed a tetrahedral mesher in FDS.

**3. Contact with 1-D solids.**
- A 1-D-backed face that abuts an FEM solid is coupled through a Robin/flux exchange equivalent to the EXPOSED backing. This removes HT3D limitation 2.

**Input sketch.**
```
&SURF ID='STEEL', HT3D=T, CELL_SIZE=0.005 /                       ! existing syntax → FEM solver
&OBST XB=..., SURF_ID='STEEL', MATL_ID='STEEL' /
&GEOM ID='BEAM', BINARY_FILE='beam.bingeom', SURF_ID='STEEL', MATL_ID='STEEL', HT3D=T /
&MATL ID='WOOD', CONDUCTIVITY_TENSOR=0.12,0.12,0.35, GRAIN_DIRECTION=0,0,1, ... /   ! new, optional
```

### B.5 Gas–solid interface

- **Faces backed by the FEM** (WALL cells of HT3D OBSTs, CFACEs of HT3D GEOMs):
  - Skip `SOLID_HEAT_TRANSFER` (the guards at `wall.f90:1875`/`1910` become a "FEM-owned" test).
  - After each FEM solve, write the Phase 0.4 outputs into `B1`. Gas-phase code is untouched.
- **Mapping.** Area-weighted and conservative in both directions:
  - Gas → facet: h, T_g, q_rad,in, ZZ_G, and q_liquid or the film node (Part C).
  - Facet → gas: TMP_F, Q_CON_F, Q_RAD_OUT, ε, mass fluxes.
  - For OBSTs the map is one-to-few by construction. For GEOM it goes through `CUT_FACE%BODTRI` and `CF%AREA`, using the Phase 0.2 machinery.
- **Cadence.**
  - Keep the `WALL_INCREMENT` cadence and `DT_BC` accumulation of `WALL_BC` (`wall.f90:89-99`).
  - Reproduce the `Q_CON_F`/`Q_RAD_OUT` time averaging (`wall.f90:2920-2928`, `2968-2969`) for energy consistency.
  - The FEM solve replaces `MESH_EXCHANGE(66)` + `HT3D_TEMPERATURE_EXCHANGE` at `main.f90:909-920` for FEM-owned solids.
- **Output.**
  - Interpolate FEM temperature to the gas-cell centres inside OBST solids for SLCF. This reuses the existing trick of writing solid temperature into the gas `TMP` (`wall.f90:3673`).
  - New point query `DEVC QUANTITY='SOLID TEMPERATURE'` (XYZ inside a solid), which generalizes `INSIDE WALL TEMPERATURE`.
  - BNDF on wall cells works unchanged; GEOM boundary data goes through BNDE.
  - Volume output for Smokeview later.

### B.6 Pyrolysis and mass transport in 3-D

- **Kinetics at nodes.** Nodal (dual-volume) component densities ρ_i. Rates come from the refactored pure kinetics routine (Phase 0.4): A, E, n_s, n_T, n_O2, `NU_*` and H_R(T).
- **Gas routing to the surface.**
  - The 1-D model releases gas instantly at its own face.
  - The 3-D equivalent routes each node's gas production to its **nearest exposed facet**, using a precomputed graph-distance field (Dijkstra from the boundary nodes) and the same instant-release assumption.
  - This replaces the `LAYER_DIVIDE`/`HT3D_WEIGHT` routing heuristics. The semantic change must be documented.
  - Darcy flow in porous char is a later option.
- **In-depth oxygen.** X_O2 at a node = X_O2 of its routed facet × exp(−d/`GAS_DIFFUSION_DEPTH`). This generalizes `wall.f90:3471-3491`.
- **Geometry change.**
  - v1: none. Density falls at fixed volume, as HT3D does today.
  - B6: **element deactivation**. When all elements inside a gas cell are consumed, remove that cell through the existing OBST `BURN_AWAY`/`REASSIGN_WALL_CELLS` machinery (`init.f90:4909`), and expose new facets.
  - GEOM burn-away stays out of scope, because cut-cell geometry is static.
- **Stay 1-D.** Liquid fuels, thermally thick particles, SPyro and TGA.

### B.7 Parallelization

- **Ownership.**
  - Elements belong to the rank that owns the gas mesh containing the element centroid.
  - Interface facets are therefore co-located with their WALL/CFACE owners, so interface data needs no extra communication.
  - Nodes shared across ranks become ghost nodes, exchanged with the Phase 0.2/0.3 machinery.
- **Assembly.** Element-by-element into HYPRE IJ, with off-process contributions handled by HYPRE at assemble time.
- **Global numbering.** Contiguous per-rank offsets, as in `GET_MATRIX_INDEXES_H`.
- **Load balance.** Geometric ownership can be unbalanced for solid-dominated cases.
  - Add optional graph partitioning (ParMETIS, or recursive coordinate bisection to avoid a new dependency), plus a facet-data redistribution step.
- **Invariance.** Results must be independent of the MPI process count to solver tolerance, and this is tested.

### B.8 Work by file

| File | Change |
|---|---|
| `sfem.f90` (new) | Mesh types (nodes, tets, materials, facets, ghost lists). Assembly, the θ-scheme and Picard. Pyrolysis driver. Gas-routing map. Output interpolation. Restart I/O. |
| `wall.f90` | `WALL_BC` skips FEM-owned faces and calls the FEM solve on the corrector cadence. Keep `PYROLYSIS` kinetics as a pure routine (Phase 0.4). Remove the HT3D sweep logic once legacy is retired. |
| `init.f90` | Build FEM regions after `INITIALIZE_WALL_ARRAY`/`FIND_WALL_BACK_INDEX` (`init.f90:27`, `3877`). Retire `INITIALIZE_HT3D_WALL_CELLS` (`init.f90:1831`) and the HT3D branches in `REALLOCATE_ONE_D_ARRAYS` (`init.f90:1621`). |
| `geom.f90` | Accept `VOLUS`/`MATL_ID` in `&GEOM`. Material tags per tetrahedron. Mesh-quality checks. Restore the constraint at `geom.f90:2322`. |
| `read.f90` | MATL `CONDUCTIVITY_TENSOR`/`GRAIN_DIRECTION`. `HT3D_SOLVER`. FEM controls (`THETA`, `PICARD_ITERATIONS`, solver tolerances). Relax the HT3D input errors that no longer apply (ERROR 437, exterior boundary; ERROR 423/424, exposed faces in each direction). |
| `main.f90` | Replace the HT3D exchange at `main.f90:909-920`. Setup after the back-wall exchange initialization (`main.f90:285-294`). |
| `imkl.f90` | Phase 0.3 wrapper. |
| `dump.f90` | Restart (`dump.f90:3681`, `3904` pattern). New solid quantities. Energy budget term for solid storage in `UPDATE_HRR`. |
| `CMakeLists.txt`, `Build/makefile` | New source files. Optional ParMETIS. |
| `Utilities/Python/fds_tetmesh` (new) | gmsh/TetGen wrapper that writes `.bingeom` with tetrahedra and material tags. |

### B.9 Phases

**B0. Refactor (1–2 person-months).** Phase 0.3 and 0.4; the kinetics routine becomes pure; bit-identical regression.

**B1. Serial linear conduction (3–4 person-months).**
- Box-to-tetrahedron meshing of OBSTs and `.bingeom` import.
- Constant properties and the Robin BC.
- *Exit:* `ht3d_slab`, `ht3d_sphere_24/48/96`, `ht3d_ibeam` (vs. the ANSYS FEM results already in the repo), and manufactured-solution convergence (second order in L2).

**B2. Nonlinear properties and radiation (2–3 person-months).**
- k(T) and c(T), multiple materials, anisotropic k.
- *Exit:* `ht3d_energy_conservation*`, `ht3d_radiation_exact`, SFPE cases 8 and 16, and the `heat_conduction_a–d` analytical cases on 1-D-like meshes.

**B3. MPI (3–4 person-months).**
- Ghost nodes, distributed HYPRE, ownership.
- *Exit:* `ht3d_2mesh_exact`, `ht3d_beam_intersection`, `ht3d_network`; partition-invariance tests; a scaling study (`ht3d_demo`).

**B4. Pyrolysis (3–5 person-months).**
- Nodal kinetics, gas routing, O₂ depth model.
- *Exit:* `ht3d_mass_conservation(_2)`; mass and energy closure; 1-D-limit agreement with `Verification/Pyrolysis/*`.

**B5. GEOM coupling (2–4 person-months).**
- CFACE ↔ triangle maps, the tetmesh utility, Smokeview output.
- *Exit:* `back_wall_test_2`-style OBST-vs-GEOM-vs-rotated-GEOM agreement. Rotation invariance is something HT3D cannot deliver.

**B6. Geometry change and retirement (3–5 person-months).**
- Element deactivation and burn-away for OBSTs, contact with 1-D solids, a P1 (diffusion) model for semi-transparent solids, and deprecation of legacy HT3D.

### B.10 Verification and validation

**Existing suite to pass with the FEM solver** (`Verification/Heat_Transfer`):
- `ht3d_*` (slab, sphere, I-beam, beam heating/intersection, 2-mesh, radiation, energy and mass conservation, pile, network, demo).
- `insulated_steel_*`, `back_wall*`, SFPE cases.

**New verification cases.**
- 3-D manufactured solutions with k(T) on unstructured tetrahedra.
- Anisotropic analytical solution, via a coordinate stretch of an isotropic solution.
- Two-material contact.
- Rotated-geometry invariance.
- MPI partition invariance.
- Exterior-boundary solids.

**Validation.**
- Sets that already use HT3D: `Validation/SNL_Walls`, `SP_AST`, `LEMTA_Spray`.
- `NIST_Composite_Beam`; `Insulation_Materials`.
- Pyrolysis sets, as a 1-D-limit regression: `Kashiwagi_Gasification`, `NIST_Polymers`, `FAA_Polymers`, `UMD_Polymers`, `Scaling_Pyrolysis`.
- The MaCFP condensed-phase cases.

### B.11 Risks and mitigations

- **Cost of global implicit solves every `WALL_INCREMENT`.**
  - Reuse the AMG setup; allow a solid time step that is a multiple of `DT_BC`; profile early (B3).
- **Mesh quality from external meshers.**
  - Quality checks, and documented gmsh settings.
- **Pyrolysis semantics change**: nearest-facet routing vs. `LAYER_DIVIDE`/`HT3D_WEIGHT`.
  - Document it, and keep legacy available during the transition.
- **Load imbalance** when solids dominate.
  - Optional repartitioning (B3).
- **Element deactivation interacting with OBST removal and pressure-matrix rebuild** (`GLOBAL_MATRIX_REASSIGN`, `main.f90:1678-1707`).
  - Defer to B6, after v1 is stable.

---

## 5. Part C: Interfacial (thin-film) solver for water on surfaces

### C.1 What exists today

**Surface droplets.** Water on surfaces is a set of ordinary Lagrangian liquid particles that stick on contact and then move with **prescribed kinematics**.
- `MOVE_ON_SOLID` (`part.f90:2583-2589`) does x = x_old + u·DT.
  - Vertical faces: `VERTICAL_VELOCITY` (default 0.5 m/s).
  - Horizontal faces: a random direction at `HORIZONTAL_VELOCITY` (default 0.2 m/s) (`read.f90:6350-6351`, `part.f90:2438-2472`).
  - GEOM faces: a down-slope direction (`part.f90:2281-2300`).
- Droplets drip off edges through geometric heuristics (`part.f90:2515-2553`).
- There is no film momentum, thickness transport, gas shear, surface tension, contact angle, splashing or boiling regime.

**Film thickness and wetted area** are cell-local estimates (`part.f90:3728-3784`). The depth is zeroed unless the droplet footprints cover the cell (`part.f90:3774`, `3781`).

**Heat transfer.**
- Each droplet exchanges heat with the gas and the wall in an implicit 3-way solve (`part.f90:4099-4144`), using flat-plate correlations with L = 1 m.
- The wall feels it only as a lagged, exponentially smoothed explicit flux, `LP_CPUA` → `Q_LIQUID_F` (`wall.f90:1893`, `1948`; boundary terms at `wall.f90:2875-2879`).

**Suppression.**
- Cooling.
- The empirical `exp(−∫E·m_w'' dt)` (`wall.f90:1105-1114`, `1340-1347`), applied only to specified-burning-rate surfaces. K never decays.

**The only Eulerian liquid-on-wall inventory** is condensate (`B1%AWM_AEROSOL`, `CONDENSATION_EVAPORATION`, `fire.f90:1922-2224`). It is WALL-only and uses a wall-function mass transfer.

**Topology.**
- There is no wall-cell adjacency anywhere.
- GEOM has triangle adjacency (`EDGE_FACES`).
- CFACEs have no CFACE ↔ CFACE or triangle → CFACE maps.

### C.2 Target model: FireFOAM feature parity

FireFOAM uses OpenFOAM's `surfaceFilmModels` (kinematic/thermo single layer) in a separate film region, coupled to the gas, to Lagrangian parcels and to the `reactingOneDim` pyrolysis region.

| FireFOAM / OpenFOAM sub-model | FDS plan | Phase |
|---|---|---|
| Kinematic single layer: depth-averaged mass and momentum, hydrostatic and capillary pressure, wall and gas shear | `film.f90` on the surface mesh | C1 |
| Thermo single layer: film energy, wall/gas heat transfer, phase change (evaporation) | film energy plus implicit film node in the solid solve | C2 |
| Film radiation (absorption, emission) | opaque film, ε ≈ 0.95, absorptance option 1−exp(−κδ) | C2 |
| Lagrangian–film interaction (`absorb`, `bounce`, `splashBai`) | absorb (C1); Bai–Gosman regime map (C4) | C1/C4 |
| Injection: `drippingInjection`, `curvatureSeparation` | drips from downward-facing facets; separation at convex edges | C4 |
| Contact-angle force and partial wetting (rivulets) | wetted-fraction model (C2); contact-angle force (C4) | C2/C4 |
| Thermocapillary force | optional | C4 |
| Film–pyrolysis coupling (mapped temperature/heat-flux BCs) | native: the film node sits in the 1-D or FEM solid's surface balance | C2/C5 |

Validation in the literature includes the Nusselt solution and FM Global vertical-film experiments (Meredith & de Vries 2013, critical heat flux for film rupture; de Vries, Meredith & Xin, panels at 7–46 kW/m²).

### C.3 Governing equations

The equations are depth-averaged over the film thickness δ, on each facet with unit normal **n** (pointing into the gas). **U** is the tangential film velocity.

**Mass.**
∂(ρ_l δ)/∂t + ∇_s·(ρ_l δ **U**) = ṁ''_imp − ṁ''_evap − ṁ''_drip/sep + ṁ''_src

**Momentum.**
∂(ρ_l δ **U**)/∂t + ∇_s·(ρ_l δ **U** **U**) = −δ ∇_s p_f + ρ_l δ **g**_t + **τ**_g − **τ**_w + **S**_imp

- **Film pressure:** p_f = −ρ_l (**g**·**n**) δ − σ ∇_s²δ.
  - On a floor, **g**·**n** = −g, so the pressure is hydrostatic and drives spreading.
  - On a ceiling the sign reverses. The Rayleigh–Taylor-unstable film then drips (C4).
- **Wall shear (laminar, half-parabolic profile):** **τ**_w = 3 μ_l **U**/δ.
- **Gas shear:** **τ**_g comes from the wall model's u_τ and the gas tangential velocity (`B2%U_TAU`, `B1%U_TANG`). `WALL_MODEL` must be called whenever a film is active; today it runs only for condensation, deposition or wall outputs (`wall.f90:163-167`).

**Energy (single film temperature T_f).**
∂(ρ_l δ h_l)/∂t + ∇_s·(ρ_l δ **U** h_l) = q''_w + q''_g + q''_rad,abs − q''_rad,emit − ṁ''_evap h_fg + ṁ''_imp h_imp − ṁ''_out h_l

- Wall-to-film conductance: h_fw ≈ 2k_l/δ.
- Gas side: the existing `HEAT_TRANSFER_COEFFICIENT`.
- Evaporation: a Spalding B-number / Sherwood formulation, reusing `GET_FILM_PROPERTIES` and `FORCED_CONVECTION_MODEL` (`func.f90:2862`, `3368`) or the wall-function mass transfer from `CONDENSATION_EVAPORATION` (`fire.f90:2153-2220`).

**Wetted fraction.**
- f_w ∈ [0, 1] from δ relative to a critical thickness. The critical thickness can be user-set, or derived from a contact angle and a minimum wetting rate.
- Fluxes split into a wet path (through the film) and a dry path (gas directly to the solid).

**Order-of-magnitude check (water, Γ = 0.05 kg/m/s on a vertical wall).**
- Nusselt thickness: δ = (3μΓ/(ρ²g))^{1/3} ≈ 0.25 mm.
- Mean velocity: ≈ 0.2 m/s.
- Wall-friction relaxation time: ρδ²/(3μ) ≈ 0.02 s.
- Film thermal time: ρc_pδ²/(2k_l) ≈ 0.2 s.
- Gravity-wave speed: √(gδ) ≈ 0.05 m/s.

**Consequences.**
- Advection is non-stiff at gas-phase time steps, so it can be explicit.
- Wall friction and the film–wall thermal exchange are comparable to or faster than DT, so they must be **implicit**.

### C.4 Discretization

- **Unknowns.** Cell-centred finite volume on the Phase 0.2 surface mesh, with facets as control volumes. Unknowns per facet: ρ_lδ, the momentum ρ_lδ**U** stored as a 3-D vector projected onto the facet tangent plane, and T_f.
- **Fluxes across shared edges.**
  - TVD upwind for advection, using FDS's existing limiters.
  - Green–Gauss gradients for p_f.
  - The capillary term ∇_s²δ uses the edge-based Laplacian.
- **Folded edges.** Where a flux crosses an edge between facets with different normals (90° Cartesian edges, GEOM creases), rotate the transported momentum about the shared edge (Rodrigues rotation).
  - Water flowing over a convex edge then continues down the next face, unless the separation criterion (C4) sheds it.
- **Positivity.** Limit edge mass outflow to the available mass. Dry facets have δ = 0 and **U** = 0.
- **Time stepping.**
  - Explicit advection, sub-cycled to a film CFL on Δx/(|**U**| + √(gδ)).
  - Wall friction and source terms are point-implicit.
- **Placement in the corrector.** After `PARTICLE_MASS_ENERGY_TRANSFER`/`MOVE_PARTICLES` and the particle exchange (`main.f90:888-905`), and before `WALL_BC` (`main.f90:911`).
- **Staircases.** On Cartesian staircase surfaces, the effective slope emerges from flow over the horizontal and vertical faces. Hydrostatic spreading drives it on the treads. This is a known limitation; GEOM gives true slopes.

### C.5 Coupling

**1. Solid: the film as an implicit node.**
- In the 1-D solver, the film adds one lumped node in front of the surface node: capacity ρ_l c_l δ, conductance f_w h_fw to the surface, with gas/radiation/evaporation fluxes on its top side.
  - It slots into the existing ghost-node Robin formulation (`RFACF`/`QDXKF`, `wall.f90:2867-2880`), extending the tridiagonal system by one row.
- In the FEM, it is one extra degree of freedom per boundary facet.
- **Cadence.** Wetted surfaces use `WALL_INCREMENT=1`, as `SURFACE_OXIDATION_MODEL` already forces (`read.f90:7310`), so the stiff exchange is resolved every gas step.
- **Energy split.**
  - The film solver transports film enthalpy between facets.
  - The solid solver owns the wall–film exchange and writes the updated T_f back.
  - Both use the same integrated flux, so energy is conserved exactly.
- **Replaces** the explicit `LP_CPUA` path for filmed surfaces. Surfaces without a film are unchanged.

**2. Gas.**
- Evaporated water enters as a surface mass flux of WATER VAPOR through the existing `M_DOT_G_PP_ADJUST` → `CALCULATE_ZZ_F` path (`wall.f90:1036-1461`). This gives the blowing velocity and `ZZ_F` consistently, exactly as pyrolysate does.
- On wetted area, the gas-side surface temperature and emissivity seen by convection and radiation (`TMP_F`, `EMISSIVITY` → `radi.f90:4220-4250`) are the f_w-weighted film and solid values.

**3. Lagrangian droplets.**
- **Impingement.** When a droplet reaches a filmed surface (the contact point in `MOVE_PARTICLES`, `part.f90:2438-2472`; `CFACE_SEARCH`, `part.f90:1996-2363`; `VENT_PARTICLE_EXTRACTION`, `part.f90:3067`), deposit PWT·m, its tangential momentum and its enthalpy into the facet, then remove the particle. The Bai–Gosman stick/rebound/spread/splash map comes in C4.
- **Release.**
  - Drips from downward-facing facets when δ exceeds a critical thickness. The drip size scales with the capillary length √(σ/ρg). Use the existing `DRIP_DIAMETER`.
  - Separation at convex edges from a force-ratio criterion.
  - Re-emit drops through the `PARTICLE_FACE_INSERT` machinery (`part.f90:603-944`).
- **Default unchanged.** Surfaces without a film keep today's surface-droplet model.

**4. Suppression.**
- On MATL-pyrolysis surfaces, suppression follows physically from cooling through the film node.
- `E_COEFFICIENT` can use film mass per area in place of the summed `LP_MPUA`, for continuity with existing inputs.

**5. Sources.**
- New `SURF FILM_MASS_FLUX` (kg/m²/s), plus a line-source option at a top edge, for wetted walls, water curtains and run-down.

**Input sketch.**
```
&FILM ID='WATER FILM', SPEC_ID='WATER VAPOR', CONTACT_ANGLE=40., CRITICAL_THICKNESS=1.E-4,
      DRIP_THICKNESS=2.E-3, SEPARATION_MODEL='FORCE_RATIO', SPLASH_MODEL='BAI_GOSMAN' /
&SURF ID='GYPSUM', MATL_ID='GYPSUM', THICKNESS=0.0127, FILM_ID='WATER FILM' /
&SURF ID='WETTED WALL', ..., FILM_ID='WATER FILM', FILM_MASS_FLUX=0.02 /
```

### C.6 MPI

- **Ownership.** A facet belongs to the rank that owns its WALL or CFACE. For GEOM triangles that straddle meshes, the owner is the mesh containing the (refined) triangle centroid.
- **Halos.** δ, **U** and T_f of neighbour facets are exchanged with the new Phase 0.2 `MESH_EXCHANGE` code. Edge fluxes are computed identically on both sides from the same states, so they stay conservative.
- **Invariance.** Results must be invariant to the MPI process count. `cascadempi.fds` is the template for a film MPI test.

### C.7 Work by file

| File | Change |
|---|---|
| `film.f90` (new) | Film state, sub-cycled transport, sources, drip/separation/splash, outputs, restart. |
| `smsh.f90` (new) | Phase 0.2 surface mesh and maps. |
| `part.f90` | Absorb-on-impact for filmed surfaces at `part.f90:2438-2472`, `1996-2363` and `3067`. Drip and separation re-emission through `PARTICLE_FACE_INSERT`. Bypass `LP_CPUA` for filmed surfaces. |
| `wall.f90` | Film node in the `SOLID_HEAT_TRANSFER` front BC (`wall.f90:2867-2880`). Evaporation mass flux into `M_DOT_G_PP_ADJUST`. `E_COEFFICIENT` from film mass (`wall.f90:1105-1114`). `WALL_MODEL` always on for filmed surfaces (`wall.f90:163-167`). |
| `main.f90` | Film step in the corrector (`main.f90:888-911`). Halo exchange. |
| `read.f90` | `&FILM` namelist. SURF `FILM_ID`/`FILM_MASS_FLUX`. Force `WALL_INCREMENT=1` when films exist. |
| `type.f90`, `func.f90` | Film types. Allocation and `PACK_*` for restart (`func.f90:4279`, `5127`). |
| `data.f90`, `dump.f90` | Outputs: `FILM THICKNESS`, `FILM VELOCITY`, `FILM TEMPERATURE`, `WETTED FRACTION`, `FILM MASS PER UNIT AREA`, `FILM EVAPORATION RATE` (BNDF/DEVC; BNDE for GEOM). Energy budget column. |

### C.8 Phases

**C1. Kinematic film on Cartesian walls, single mesh (4–5 person-months).**
- Surface mesh (Phase 0.2 Cartesian builder), mass and momentum, and droplet absorption.
- *Exit:* Nusselt film on inclined and vertical plates; Huppert's (1982) viscous gravity-current similarity solution on a horizontal plate; mass conservation around a box (over the top, down the sides, off the bottom).

**C2. Thermal film (3–4 person-months).**
- Energy, evaporation, the implicit film node in the 1-D solver, radiation, wetted fraction.
- *Exit:* analytical film heating on an isothermal wall; evaporation-rate check against the correlation; comparison of `water_evaporation_6/7` and `LEMTA_Spray_Cooling_*` against the surface-droplet model and the data.

**C3. MPI and GEOM (3–4 person-months).**
- Halos; GEOM triangle film mesh with refinement; CFACE maps.
- *Exit:* MPI invariance; inclined GEOM plate equals a rotated OBST staircase in the limit; `geom_particle_cascade*` with a film.

**C4. Interfacial sub-models (3–5 person-months).**
- Dripping (validate the Rayleigh–Taylor drip spacing ≈ 2π√2·l_c under a ceiling).
- Edge separation; Bai–Gosman splashing; contact-angle and partial wetting (rivulets); nucleate/film boiling curve and Leidenfrost cut-off.

**C5. Suppression and FEM surfaces (2–3 person-months).**
- Film on pyrolysing MATL surfaces and on FEM solids (after B5).
- *Exit:* `Validation/UL_NFPRF` without a tuned `E_COEFFICIENT`; FM vertical-film rupture and panel data (Meredith & de Vries; de Vries, Meredith & Xin); `Validation/LEMTA_Spray` (HT3D steel plate, so it exercises FEM + film together).

### C.9 Risks and mitigations

- **Staircase OBST geometry gives non-physical slope flow.**
  - Document it; recommend GEOM for sloped wetted surfaces.
- **CFACE fragments and merged cut cells complicate the facet ↔ CFACE map.**
  - Run the film on GEOM triangles, not CFACEs, and use area-weighted many-to-many maps.
- **Partial-wetting physics is empirical.**
  - Start with a calibrated critical-thickness model, then add the contact-angle force; validate against FM data.
- **Stiff thermal coupling.**
  - Use the implicit film node with `WALL_INCREMENT=1`. Explicit coupling would oscillate at the ≈0.2 s film time scale.
- **Interaction with OBST create/remove and burn-away.**
  - v1 forbids films on removable OBSTs. Later, rebuild the surface mesh in `REASSIGN_WALL_CELLS`.

---

## 6. Integrated roadmap

```
Phase 0 ── 0.1 MPI_COMM_FDS ───────────────┬──────────────► A2 (online WRF)
        ├─ 0.3 sparse solver wrapper ──────┼──► B1 ─► B2 ─► B3 ─► B4 ─► B5 ─► B6
        ├─ 0.4 exchange contract, kinetics ┼──► B1, C2                 │
        └─ 0.2 surface-mesh layer ─────────┴──► C1 ─► C2 ─► C3 ─► C4 ─► C5 ◄┘ (film on FEM surfaces)
A1 (offline WRF → FDS) can start immediately ─► A2a ─► A2b ─► A3
```

The schedule below assumes three developers, one per track, and shares Phase 0 between them. It runs to roughly 18 months.

| Quarter | Track A (WRF) | Track B (FEM) | Track C (film) |
|---|---|---|---|
| Q1 | Phase 0.1; A1 preprocessor and boundary arrays | Phase 0.3/0.4; B1 | Phase 0.2 (Cartesian builder) |
| Q2 | A1 complete; Askervein / WRF-LES ABL | B2 | C1 |
| Q3 | A2a (MPMD, SFIRE-equivalent mode) | B3 | C2 |
| Q4 | A2b (LES nest, `SOURCES` feedback) | B4 | C3 |
| Q5 | A3 (FireFlux, sensitivity) | B5 | C4 |
| Q6 | `STATE` feedback, fuel moisture, docs | B6, legacy HT3D deprecation | C5, docs |

---

## 7. Decisions needed before starting

1. **WRF flavour and version.** NCAR WRF v4.x (WRF-Fire in `phys/`) or OpenWFM WRF-SFIRE. This decides where the WRF-side hooks go and which flux arrays are reused.
2. **Coupler backend.** OASIS3-MCT (recommended if the target WRF has its OASIS interface) or a bespoke MPI intercommunicator.
3. **Scope of "replace".** Replace HT3D and keep the 1-D solver (recommended), or also replace the 1-D solver. Replacing it as well is not recommended (Section B.2).
4. **Element types.** P1 tetrahedra only (recommended for v1), or also Q1 hexahedra for OBSTs.
5. **Tooling dependencies.** Python preprocessors in `Utilities/` (`wrf2fds` with netCDF/pyproj; `fds_tetmesh` with gmsh) keep the FDS runtime dependency-free beyond HYPRE and an optional OASIS3-MCT. Is that acceptable?
6. **Upstream intent.** If the goal is to merge into `firemodels/fds`, engage the NIST developers early. This matters for the 2014 `MEAN_FORCING`/WRF history, HT3D ownership, and manual and Firebot conventions.
7. **Priority** if resources are limited. Each track is independent after Phase 0.

---

## 8. Pre-existing issues found during the survey

These were found while mapping the code. They are not part of the three features, but they affect them.

1. **Droplet wall cooling is undercounted when the evaporation step sub-cycles.**
   - `LP_CPUA` is accumulated from `Q_CON_WALL/DT` (`part.f90:4361-4362`). `Q_CON_WALL` is only the *last* sub-step's energy.
   - `Q_CON_SUM` is accumulated over all sub-steps (`part.f90:4338`) but never used.
   - The gas-phase energy budget (`Q_DOT(5)`, `part.f90:4322`) counts every sub-step. So the wall cooling and the budget disagree in hot-wall cases.
2. **Documented minimum film thickness is not implemented.** The Technical Reference Guide's droplet film thickness, δ = max(δ_min = 10⁻⁵ m, …) (`Particle_Chapter.tex`, "Droplet Transport on a Surface"), does not match the code. The code uses a footprint-coverage test instead (`part.f90:3774`, `3781`).
3. **Outdated Verification Guide note.** It says HT3D "is not yet connected to the 1-D pyrolysis" (`FDS_Verification_Guide.tex:4868`). `Verification/Heat_Transfer/ht3d_mass_conservation.fds` exercises HT3D with pyrolysis.
4. **`LEVEL_SET_MODE=5` is undocumented, and its prescribed-wind branch is dead.**
   - Mode 5 sets `LEVEL_SET_COUPLED_WIND=.TRUE.` (`read.f90:2012`).
   - Its prescribed-wind branch (`vege.f90:356-365`) is then overwritten by the CFD sampling block (`vege.f90:371-402`).
   - It is relevant to A2a, which needs a fire-coupled mode with externally supplied spread-rate wind.
5. **The lighter back-wall exchange is never used.** `MESH_EXCHANGE(67)` (`PACK_WALL_LITE`) is never called. Code 66 sends the full `PACK_WALL` every step (`main.f90:1009`). This is a performance item relevant to B3 and C3.
6. **HYPRE solver settings are compile-time parameters** (`imkl.f90:204-229`). Phase 0.3 addresses this.

---

## 9. References

- Mandel, J., Beezley, J. D., & Kochanski, A. K. (2011). Coupled atmosphere-wildland fire modeling with WRF 3.3 and SFIRE 2011. *Geosci. Model Dev.*, 4, 591–610. <https://gmd.copernicus.org/articles/4/591/2011/>
- Kochanski, A. K. et al. Evaluation of WRF-SFIRE performance with field observations from the FireFlux experiment. <https://arxiv.org/abs/1206.3345>
- Lin, D. et al. (2021). WRF4PALM v1.0: a mesoscale dynamical driver for the microscale PALM model system 6.0. *Geosci. Model Dev.*, 14, 2503–2524. <https://gmd.copernicus.org/articles/14/2503/2021/>
- Kadasch, E. et al. (2021). Mesoscale nesting interface of the PALM model system 6.0. *Geosci. Model Dev.*, 14, 5435–5465. <https://gmd.copernicus.org/articles/14/5435/2021/>
- Muñoz-Esparza, D., Kosović, B., et al. Stochastic cell perturbation method for mesoscale-to-microscale (WRF → LES) inflow turbulence (Boundary-Layer Meteorology, 2014; JAMES, 2015).
- OpenFOAM `surfaceFilmModels` (`kinematicSingleLayer`, `thermoSingleLayer`). <https://api.openfoam.com/2312/classFoam_1_1regionModels_1_1surfaceFilmModels_1_1thermoSingleLayer.html>
- FireFOAM overview. <https://openfoamwiki.net/index.php/FireFoam>
- Meredith, K. V., & de Vries, J. FireFOAM thin-film model and vertical-film validation (OpenFOAM Workshop 2010; critical heat flux for film rupture, 2013). <https://www.tfd.chalmers.se/~hani/OFW5/timeSchedule_workshop_online-filer/Abstract_final/KarlMeredithAbstractOFW5.pdf>
- Bai, C., & Gosman, A. D. (1995). Development of methodology for spray impingement simulation. SAE 950283.
- Huppert, H. E. (1982). The propagation of two-dimensional and axisymmetric viscous gravity currents over a rigid horizontal surface. *J. Fluid Mech.*, 121, 43–58.
