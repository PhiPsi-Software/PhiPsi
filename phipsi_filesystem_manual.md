# PhiPsi File System Manual
+ This manual describes the file formats (input files, output files, restart files) that the PhiPsi solver reads and writes, plus the file-to-visualization mapping used by the companion PPView post-processor.
+ Author: Fang Shi
+ Website: http://phipsi.top
+ Email: shifang@hyit.edu.cn / shifang@ustc.edu.cn
+ Updated on 2026-06-29.

---

Click [here](http://phipsi.top/source/PhiPsi_Win64_Latest.zip) to download the latest version of PhiPsi.

Click [here](https://sourceforge.net/projects/phipsi/files/) to download the latest version of PPView.

For keyword (`.kpp`) reference, see the companion [PhiPsi Keywords Manual](http://www.phipsi.top/phipsi_keywords_manual.html).

---

# Contents of this file system manual

[File naming rules](#section-File-naming-rules)

# Input files

[Keywords input (.kpp)](#section-Keywords-input-kpp)

[Required mesh files (.node, .elem)](#section-Required-mesh-files-node-elem)

[Boundary-condition files (.boux/.bouy/.bouz, .buxn/.buyn/.buzn)](#section-Boundary-condition-files-bouxx-buyn-buzn)

[Force / point-load files (.focx/.focy/.focz)](#section-Force-point-load-files-focx-focy-focz)

[Initial-condition files for dynamics (.ive*, .iac*, .idp*)](#section-Initial-condition-files-ive-iac-idp)

[Degree-of-freedom coupling files (.dofx/.dofy)](#section-DOF-coupling-files-dofx-dofy)

[Field-problem input files (.fbvl, .fbqn, .fbiv)](#section-Field-problem-files-fbvl-fbqn-fbiv)

[Hydraulic-fracturing wellbore files (.bhpc)](#section-HF-wellbore-files-bhpc)

[Earthquake node list (.eqnl)](#section-Earthquake-node-list-eqnl)

[3D model boundary-detection outputs (.outl/.outa/.outn)](#section-3D-boundary-detection-outputs-outl-outa-outn)

[PPView geometry metadata (.boud, .focd)](#section-PPView-geometry-metadata-boud-focd)

# Output files

[VTK files (.vtk, _CRACK_.vtk)](#section-VTK-files)

[Displacement files (.disp_, .disn_, .dien_, .disg_, .dipc_, .dinc_, .veln_, .acln_)](#section-Displacement-files)

[Stress files (.strn_, .strg_, .stnc_, .sttn_, .elss_)](#section-Stress-files)

[Strain files (.sran_, .srac_)](#section-Strain-files)

[Gauss-point files (.gcor_, .disg_, .damg_, .elgn_, .strg_)](#section-Gauss-point-files)

[SIF / fracture-energy files (.sifs_, .sift_, .prst_, .enee)](#section-SIF-fracture-energy-files)

[Crack coordinates (.crax/.cray/.craz, .crxo/.cryo, .cnox/.cnoy/.cnoz)](#section-Crack-coordinates)

[Crack meshes (.cms1/.cms2/.cms3, .cmso, .cmap, .cape, .capf, .ctap)](#section-Crack-meshes)

[Crack vertex and tangent files (.cvpx/.cvpy/.cvpz, .cvxx/.cvxy/.cvxz, ...)](#section-Crack-vertex-tangent-files)

[Crack calculation-point files (.apex/.apey, .cori, .cpre, .cvel, .cqua, .cohx/.cohy, .crrd)](#section-Crack-calculation-point-files)

[XFEM enrichment files (.ennd, .elty, .posi, .elts, .posh, .ennh, ...)](#section-XFEM-enrichment-files)

[Crack-tip and baseline files (.celt, .celv, .celj, .celc, .ctty, .blab, .blvx/.blvy/.blvz, .tere)](#section-Crack-tip-baseline-files)

[Crack-node S1 vector files (.cndx/.cndy/.cndz)](#section-Crack-node-S1-vector-files)

[Hole and inclusion files (.hlcr, .ehcr, .jzcr, .jzpx/.jzpy)](#section-Hole-inclusion-files)

[Cross-interface files (.cscr)](#section-Cross-interface-files-cscr)

[Natural-fracture files (.nfcx/.nfcy/.nfcz, .ncrx/.ncry)](#section-Natural-fracture-files)

[Energy history files (.ener, .edye, .edtm, .edap)](#section-Energy-history-files)

[Crack length history (.dcrl, .dprl, .idtm, .iite, .ilth, .ipre)](#section-Crack-length-history)

[Dynamic files (.idtm, .edtm, .dcrl, .dprl)](#section-Dynamic-files)

[Lumped mass files (.lpmf_, .lpmx_, .lpma_)](#section-Lumped-mass-files)

[Hydraulic-fracturing output files (.hftm, .injp, .wbfp, .wbpt, .ihft, .icpt, .icpt_LS, .icpt_CT)](#section-HF-output-files)

[Proppant output files (.ccon_, .pokf_, .cond_, .wpnp_, .epcr_, .Saved_Filename)](#section-Proppant-output-files)

[Contact, cohesive and element-state files (.elcs_, .elco_, .elss_, .kiel_)](#section-Contact-cohesive-element-state)

[DOF / surface-load force files (.fxdf_, .fydf_, .fzdf_, .fxsl_, .fysl_, .fzsl_)](#section-DOF-surface-load-files)

[Stiffness-matrix files (.skxf, .csrn_, .csra_, .csrj_, .csri_)](#section-Stiffness-matrix-files)

[3D fluid-crack files (.fenx/.feny/.fenz, .fnnx/.fnny/.fnnz, .fnux/.fnuy/.fnuz, .fnlx/.fnly/.fnlz, .fexx/.fexy/.fexz/.feyx/.feyy/.feyz/.fezx/.fezy/.fezz, .cpfn_, .cpno_, .ccpx/.ccpy/.ccpz, .cnlx/.cnly/.cnlz, .cmse_)](#section-3D-fluid-crack-files)

[Miscellaneous files (.fraz, .seed, .post, .fdcu, .fccu, .rbco_, .fdvl, .ecfv, .sccx/.sccy/.sccz, .scdx/.scdy/.scdz)](#section-Miscellaneous-files)

[Console / log files (PhiPsi_Console_Window.log, current_folder.dat)](#section-Console-log-files)

# PPView file-reading map

[Cloud plot / contour rendering](#section-Cloud-plot-rendering)

[Curve plot / line chart rendering](#section-Curve-plot-rendering)

[Vector plot / arrow rendering](#section-Vector-plot-rendering)

[Crack visualization (2D and 3D)](#section-Crack-visualization)

[Animation playback](#section-Animation-playback)

# Data format examples

[Example A — 2D simple FEM (2d_fem)](#section-Example-A-2d-fem)

[Example B — 2D crack propagation (2d_L_shaped_specimen)](#section-Example-B-2d-L-shaped-specimen)

[Example C — 3D uniaxial tension (3d_uniaxial_tension)](#section-Example-C-3d-uniaxial-tension)

[Example D — 3D hydraulic fracturing (3d_hf_toughness_dominated)](#section-Example-D-3d-hf-toughness-dominated)

# Quick-reference summary tables

[All input-file summary table](#section-All-input-files-summary)

[All output-file summary table](#section-All-output-files-summary)

---

<a id="section-File-naming-rules"></a>

# File naming rules

+ All PhiPsi input and output files share a common base path, which is built by joining the two keywords `*Work_Directory` and `*Filename` from the `.kpp` input. Every file in this manual has the form `<basename>.{ext}[_{suffix}]`.

 ~~~bash
 *Work_Directory
 X:\PhiPsi_Work\FEM

 *Filename
 FEM
 ~~~

  produces `basename = X:\PhiPsi_Work\FEM\FEM`. From that base, the solver reads `FEM.node`, `FEM.elem`, `FEM.boux`, … and writes `FEM.disp_1`, `FEM.disn_1`, `FEM.sifs_1`, `FEM.ener`, ….

+ Input files (read by the solver) use a fixed extension pattern; they have **no** step-number suffix. The extension alone selects what kind of data is in the file.

+ Output files (written by the solver) usually carry a **per-step suffix** `_i` where `i` is the load-step (substep) index. For VTK files the step index is zero-padded to 5 digits (`_00001.vtk`); for other outputs the step number is written without leading zeros (`_1`, `_2`, …). A few output files do **not** carry a step suffix — they are written once (e.g. `.ener`, `.seed`, `.post`, `.wbfp`) or once per analysis type (e.g. `.hlcr` for holes).

+ Some output files can be written as plain ASCII text or as a binary stream. The choice is governed by `*Key_Data_Format` (1 = ASCII, 2 = binary). When a file is documented as "Format: ASCII / Binary", the binary stream mirrors the ASCII column layout but with no fixed-width formatting — values are written as raw floating-point records.

---

# Input files

<a id="section-Keywords-input-kpp"></a>

## Keywords input (.kpp)

+ ***<font color=green>.kpp</font>** - Keywords input file. Drives every PhiPsi run; selects analysis type, mesh files, material, solver, sub-steps, and which outputs are produced.

  - **Pattern**: `<Work_Directory>\*Work_Directory points here*.kpp`
  - **Required**: Yes — solver aborts if either `*Work_Directory` or `*Filename` is missing.
  - **Format**: Plain text, `*Keyword` lines followed by their value(s); lines starting with `%` are comments.
  - **Dimensions**: both
  - **Convention**: Keywords are case-insensitive; supports parameter definitions (`num_of_crack = 2`) and basic arithmetic (`1+1`, `(15e6-10e6)*(3-1)`).
  - **Reference**: See [PhiPsi Keywords Manual](http://www.phipsi.top/phipsi_keywords_manual.html) for the complete list of 200+ keywords.

  Example:
  ~~~bash
  % Working directory and filename.
  *Work_Directory
  X:\PhiPsi_Project\Test_Examples\2d_fem

  *Filename
  fem_2d

  *Key_Analysis_Type
  1                 % Quasi-static.

  *Key_Type_2D
  1                 % Plane stress.

  *Key_SLOE
  6                 % MUMPS solver.

  *Num_Substeps
  1

  *Material_Para_1
  70.0e9,0.3,2700.0,1.0,1.0e6,1.0e6,100.0e6,0.0,...
  ~~~

<a id="section-Required-mesh-files-node-elem"></a>

## Required mesh files (.node, .elem)

+ ***<font color=green>.node</font>** - Node coordinate file.

  - **Pattern**: `<basename>.node`
  - **Required**: Yes — solver aborts if missing.
  - **Dimensions**: 2D = `X Y`; 3D = `X Y Z`.
  - **Format**: ASCII, whitespace-separated real numbers; one node per line; no header. Read by `Tool_Read_File` after `Tool_Count_Lines` counts rows.
  - **Units**: meters (SI) or millimeters (mm-ton-s) depending on `*Key_Unit_System`.

  Example (`Test_Examples/2d_L_shaped_specimen/2d_L_shaped_specimen.node`):
  ~~~bash
    0.00000000E+00   0.00000000E+00 
    0.25000000E+00   0.00000000E+00 
    0.16666667E-01   0.00000000E+00 
    0.33333333E-01   0.00000000E+00 
    0.50000000E-01   0.00000000E+00 
    0.66666667E-01   0.00000000E+00 
    0.83333333E-01   0.00000000E+00 
    0.10000000E+00   0.00000000E+00 
    0.11666667E+00   0.00000000E+00 
    0.13333333E+00   0.00000000E+00 
  ~~~

  3D example (`Test_Examples/3d_uniaxial_tension/3d_uniaxial_tension.node`):
  ~~~bash
    0.00000000E+00   0.00000000E+00   0.00000000E+00
    0.11000000E+02   0.00000000E+00   0.00000000E+00
    0.36666667E+00   0.00000000E+00   0.00000000E+00
    0.73333333E+00   0.00000000E+00   0.00000000E+00
    0.11000000E+01   0.00000000E+00   0.00000000E+00
    0.14666667E+01   0.00000000E+00   0.00000000E+00
    0.18333333E+01   0.00000000E+00   0.00000000E+00
    0.22000000E+01   0.00000000E+00   0.00000000E+00
  ~~~

+ ***<font color=green>.elem</font>** - Element connectivity and material-number file.

  - **Pattern**: `<basename>.elem`
  - **Required**: Yes.
  - **Dimensions**: 2D = 5 columns (`N1 N2 N3 N4 Mat`); 3D = 9 columns (`N1 N2 N3 N4 N5 N6 N7 N8 Mat`).
  - **Format**: ASCII, whitespace-separated integers; one element per line; no header.
  - **Units**: node indices are 1-based integers; `Mat` references a `*Material_Para_<n>` keyword in the `.kpp`.
  - **Element type**: 2D = 4-node quad (linear); 3D = 8-node hexahedron (linear). Triangular 6-node elements are also supported by some analyses (see VTK cell types in the VTK output section).

  2D example (`Test_Examples/2d_fem/fem_2d.elem`):
  ~~~bash
        1.       3.      81.      62.       1.
        3.       4.     100.      81.       1.
        4.       5.     119.     100.       1.
        5.       6.     138.     119.       1.
        6.       7.     157.     138.       1.
        7.       8.     176.     157.       1.
        8.       9.     195.     176.       1.
        9.      10.     214.     195.       1.
       10.      11.     233.     214.       1.
       11.      12.     252.     233.       1.
  ~~~

  3D example (`Test_Examples/3d_uniaxial_tension/3d_uniaxial_tension.elem`):
  ~~~bash
        1.       3.     121.     120.    1072.    1073.    6723.    4722.       1.
        3.       4.     150.     121.    1073.    1113.    7883.    6723.       1.
        4.       5.     179.     150.    1113.    1153.    9043.    7883.       1.
        5.       6.     208.     179.    1153.    1193.   10203.    9043.       1.
        6.       7.     237.     208.    1193.    1233.   11363.   10203.       1.
        7.       8.     266.     237.    1233.    1273.   12523.   11363.       1.
  ~~~

<a id="section-Boundary-condition-files-bouxx-buyn-buzn"></a>

## Boundary-condition files (.boux/.bouy/.bouz, .buxn/.buyn/.buzn)

These six files declare nodal Dirichlet boundary conditions.

+ ***<font color=green>.boux</font>** - X-direction zero-displacement (fixed DOF) node list.

  - **Pattern**: `<basename>.boux`
  - **Required**: No (warned if missing on a problem with intended fixed BCs).
  - **Dimensions**: 2D and 3D.
  - **Format**: ASCII, one node index per line; no header.
  - **Reading**: `Num_Bou_x = Tool_Count_Lines(temp_name)`; then `Bou_x = int(Temp_DATA(:,1))`.

  Example (`Test_Examples/2d_fem/fem_2d.boux`):
  ~~~bash
        1.
       42.
       62.
       63.
       64.
       65.
       66.
       67.
       68.
       69.
  ~~~

+ ***<font color=green>.bouy</font>** - Y-direction zero-displacement node list. Same format as `.boux`. In 2D and 3D.

  Example (`Test_Examples/2d_fem/fem_2d.bouy`):
  ~~~bash
        1.
       42.
       62.
       63.
       64.
       65.
       66.
       67.
       68.
       69.
  ~~~

+ ***<font color=green>.bouz</font>** - Z-direction zero-displacement node list. **3D only**.

  - **Format**: identical to `.boux`/`.bouy`.

+ ***<font color=green>.buxn</font>** - Non-zero prescribed X displacement. Format: `Node_ID  Value` (2 columns).

+ ***<font color=green>.buyn</font>** - Non-zero prescribed Y displacement. Same as `.buxn`.

+ ***<font color=green>.buzn</font>** - Non-zero prescribed Z displacement. 3D only.

  - **Required**: No; if absent, only zero BCs are read.
  - **Format**: ASCII, `node_id value` per line.
  - **Units**: same as `.node` (m or mm).

<a id="section-Force-point-load-files-focx-focy-focz"></a>

## Force / point-load files (.focx/.focy/.focz)

+ ***<font color=green>.focx</font>** - Nodal point load in X direction.

  - **Pattern**: `<basename>.focx`
  - **Required**: No.
  - **Format**: ASCII, `Node_ID  Value` per line (2 columns); no header.
  - **Units**: Newtons in SI; tonnes of force per node in mm-ton-s.

  Example (`Test_Examples/2d_L_shaped_specimen/2d_L_shaped_specimen.focx` — single row):
  ~~~bash
    0.137028414927E-04  0.247500000000E-01
  ~~~

+ ***<font color=green>.focy</font>** - Y-direction nodal force. Same format.

  Example (`Test_Examples/2d_fem/fem_2d.focy`):
  ~~~bash
       22.   0.22500000E+05
      462.   0.22500000E+05
      482.   0.45000000E+05
      483.   0.45000000E+05
      484.   0.45000000E+05
      485.   0.45000000E+05
      486.   0.45000000E+05
      487.   0.45000000E+05
      488.   0.45000000E+05
      489.   0.45000000E+05
  ~~~

+ ***<font color=green>.focz</font>** - Z-direction nodal force. **3D only**. Same 2-column format.

  3D example (`Test_Examples/3d_uniaxial_tension/3d_uniaxial_tension.focz`):
  ~~~bash
      962.   0.33586629E+06
     1003.   0.33586629E+06
     1004.   0.67173257E+06
     1005.   0.67173257E+06
     1006.   0.67173257E+06
  ~~~

<a id="section-Initial-condition-files-ive-iac-idp"></a>

## Initial-condition files for dynamics (.ive*, .iac*, .idp*)

These nine files supply initial nodal velocity / acceleration / displacement for implicit-dynamic (`*Key_Analysis_Type = 2`) and explicit-dynamic analyses.

+ ***<font color=green>.ivex, .ivey, .ivez</font>** - Initial nodal velocity per DOF.
+ ***<font color=green>.iacx, .iacy, .iacz</font>** - Initial nodal acceleration per DOF.
+ ***<font color=green>.idpx, .idpy, .idpz</font>** - Initial nodal displacement per DOF (3D only via `Read_Geo_3D.F90`; for 2D use `.buxn`/`.buyn` instead).

  - **Pattern**: `<basename>.{ive|iac|idp}{x|y|z}`
  - **Required**: No.
  - **Format**: ASCII, `Node_ID  Value` per line; 2 columns.
  - **Units**: m/s for velocity, m/s² for acceleration, m for displacement.

<a id="section-DOF-coupling-files-dofx-dofy"></a>

## Degree-of-freedom coupling files (.dofx/.dofy)

+ ***<font color=green>.dofx</font>** - Multi-point-constraint coupling set in X DOF.
+ ***<font color=green>.dofy</font>** - Multi-point-constraint coupling set in Y DOF.

  - **Format**: ASCII, `Set_Number  Node_ID` per line (2 columns). Nodes with the same `Set_Number` are tied together in the X (or Y) DOF.
  - **Required**: No.
  - **Example**: `1 5\n1 6\n2 10` ties nodes 5 and 6 together in the coupled group 1; node 10 is alone in group 2.

<a id="section-Field-problem-files-fbvl-fbqn-fbiv"></a>

## Field-problem input files (.fbvl, .fbqn, .fbiv)

Used with `*Key_Analysis_Type = 15` (field problem such as heat transfer or pore pressure). 2D-only.

+ ***<font color=green>.fbvl</font>** - Field-problem boundary value list. Each row: `Node_ID  Value`. Zero values are treated as fixed-Dirichlet; non-zero values are applied as a known scalar (e.g. prescribed temperature).
+ ***<font color=green>.fbqn</font>** - Boundary flux (Neumann). Same 2-column format. Units depend on the field variable (W/m² for heat, m/s for pore-flow).
+ ***<font color=green>.fbiv</font>** - Initial value per node. Same 2-column format. Sets the starting scalar value at each listed node.

<a id="section-HF-wellbore-files-bhpc"></a>

## Hydraulic-fracturing wellbore files (.bhpc)

+ ***<font color=green>.bhpc</font>** - Bottom-hole-pressure time curve.

  - **Pattern**: `<basename>.bhpc`
  - **Required**: Conditional — read only when `*Key_Analysis_Type` is 16 or 17 (wellbore analysis), `*Key_Gas_Production = 1`, and `*Key_Changing_BHP = 1`.
  - **Format**: ASCII, two columns `Time  Pressure` per row.
  - **Units**: seconds and Pa.

<a id="section-Earthquake-node-list-eqnl"></a>

## Earthquake node list (.eqnl)

+ ***<font color=green>.eqnl</font>** - List of node IDs to apply earthquake acceleration.

  - **Pattern**: `<basename>.eqnl`
  - **Required**: Conditional — read only when `*EQ_Ac_nodes_list_method = 2` (list method).
  - **Format**: ASCII, one node index per line; no header.
  - **Limit**: up to 5000 nodes.

<a id="section-3D-boundary-detection-outputs-outl-outa-outn"></a>

## 3D model boundary-detection outputs (.outl/.outa/.outn)

These are **written by** the solver (in `Read_Geo_3D.F90` and helpers) as a side-effect of importing the 3D geometry; they are not user input but appear in the working directory after a 3D run.

+ ***<font color=green>.outl</font>** - 3D model outer boundary edge list. Format: `Node1  Node2` per row (integers).
+ ***<font color=green>.outa</font>** - 3D model outer surface-area edges. Same 2-column format.
+ ***<font color=green>.outn</font>** - 3D model exterior surface node numbers. Format: one node index per line.

<a id="section-PPView-geometry-metadata-boud-focd"></a>

## PPView geometry metadata (.boud, .focd)

These pipe-delimited text files are written by PPView when the user builds the model graphically; they mirror what `.boux`/`.bouy`/`.focx`/`.focy` mean but can express boundary-conditions and loads as functions of coordinate ranges or single coordinates.

+ ***<font color=green>.boud</font>** - Boundary-condition metadata.

  - **Format**: one pipe-delimited record per line:
    ~~~bash
    fix_nodes|X|1,2,3,4,5
    fix_coord_single|0|0.0|0.001|All
    fix_coords_range|0|10|0|5|0|0|0.001|X|Y
    fix_all_boundaries|left
    ~~~

+ ***<font color=green>.focd</font>** - Applied-force metadata.

  - **Format**: pipe-delimited:
    ~~~bash
    force_nodes|Y|1000.0|22,462,482
    force_coord_single|y|0.05|0.001|0|225000|0
    pressure_timevar|0|10|0|5|0|0|0.001|csv_file.csv
    init_vel|0|10|0|5|0|0|0.001|0|0|-9.8
    ~~~

  PPView translates these records back into the simple `.boux`/`.bouy`/`.focx`/`.focy`/`.focz` ASCII lists before invoking the solver.

---

# Output files

<a id="section-VTK-files"></a>

## VTK files (.vtk, _CRACK_.vtk)

The legacy ASCII VTK files written by `Save_vtk_file.f90` and `Save_vtk_file_for_Crack.f90`. Both are enabled by `*Key_Save_vtk = 1`.

+ ***<font color=green>{Full_Pathname}_NNNNN.vtk</font>** - Main solution snapshot at substep `NNNNN`.

  - **Pattern**: `<basename>_<isub5>.vtk` where `<isub5>` is the 5-digit zero-padded substep index (e.g. `00001`).
  - **Required**: No — gated by `*Key_Save_vtk`.
  - **Format**: Legacy ASCII VTK (`# vtk DataFile Version 4.0` header).
  - **Contents** (in order):
    - `POINTS Num_Node double` — `(X, Y, Z)` triples, `F12.6` per coordinate. For 2D, Z is `0.0`.
    - `CELLS Num_Elem ind` — element connectivity, **0-indexed** (1 is subtracted from each node ID before writing). For 2D quad: `4 n1 n2 n3 n4`. For 3D hex: `8 n1 n2 … n8`.
    - `CELL_TYPES Num_Elem` — VTK cell-type integer per element (5 = tri, 9 = quad, 10 = tetrahedron, 12 = hexahedron, 22 = quadratic tri, etc.).
    - `CELL_DATA Num_Elem` block:
      - `SCALARS Material_ID int` and `LOOKUP_TABLE default`
      - `SCALARS Element_ID int` and `LOOKUP_TABLE default`
      - `SCALARS Element_Permeability_{xx|yy|zz|xy|yz|xz} double` (if pore-pressure elements allocated)
      - `SCALARS Enriched_Element_Type int` (if `Enrich_Freedom > 0`)
    - `POINT_DATA Num_Node` block:
      - `VECTORS Displacement double` — `(Ux, Uy, Uz)` per node, `E20.12`.
      - `SCALARS stress_xx|yy|xy|vm double` for 2D; `stress_xx|yy|zz|xy|yz|xz` for 3D (if allocated).
      - `SCALARS Node_Number integer`.
      - `SCALARS Enriched_Node_Type int` (if enriched).

  Example (`Test_Examples/2d_fem/fem_2d_00001.vtk`):
  ~~~bash
  # vtk DataFile Version 4.0
  X:\PhiPsi_Project\Test_Examples\2d_fem\fem_2d  - results from increment 00001
  ASCII
  DATASET UNSTRUCTURED_GRID

  POINTS 861 double
      0.050000    0.000000    0.000000
     -0.040000    0.030000    0.000000
      0.049610    0.006229    0.000000
      0.048448    0.012361    0.000000
      0.046531    0.018300    0.000000
      0.043888    0.023955    0.000000
      0.040562    0.029236    0.000000
      0.036604    0.034061    0.000000
      0.032075    0.038356    0.000000
      0.027047    0.042053    0.000000
      0.021598    0.045095    0.000000
      0.015811    0.047434    0.000000
      0.009779    0.049034    0.000000
      0.003594    0.049871    0.000000
     -0.002647    0.049930    0.000000
  ~~~

  - **Compatibility**: Open in ParaView, VisIt, F3D, or any VTK-compatible viewer. The cell-connectivity numbering is **0-based** (1 is subtracted from each node ID before writing).

+ ***<font color=green>{Full_Pathname}_CRACK_NNNNN.vtk</font>** - Crack geometry only (separate file for visualisation clarity).

  - **Pattern**: `<basename>_CRACK_<isub5>.vtk`.
  - **Contents**:
    - 2D: `VTK_LINE` cells (cell type 3) — `2 n1 n2` per row.
    - 3D: `VTK_TRIANGLE` cells (cell type 5) — `3 n1 n2 n3` per row.
    - `SCALARS Crack_Node_Aperture double` (3D only) per crack node.
    - `SCALARS Crack_Node_Number integer`.
    - `SCALARS Crack_Element_ID int`.

<a id="section-Displacement-files"></a>

## Displacement files (.disp_, .disn_, .dien_, .disg_, .dipc_, .dinc_, .veln_, .acln_)

All written by `Save_Disp.f90`, `Save_Velocity_Accel.f90`, and the related routines. ASCII vs binary is controlled by `*Key_Data_Format`.

+ ***<font color=green>.disp_<i></font>** - Global displacement vector in raw global-DOF order.

  - **Pattern**: `<basename>.disp_<isub>`.
  - **Format**: ASCII = one `(E20.12)` value per line; binary = unformatted stream of floating-point records (8 bytes per value for double precision).
  - **Columns**: 1 (raw `DISP(i)` for i = 1 … Total_FD).
  - **Ordering**: `Ux1, Uy1, Uz1, Ux2, Uy2, Uz2, …` for 3D; `Ux1, Uy1, Ux2, Uy2, …` for 2D.
  - **Units**: SI (m). Not unit-system-scaled at write time.

+ ***<font color=green>.disn_<i></font>** - Nodal displacements indexed by node number — the user-friendly companion to `.disp_`.

  - **Pattern**: `<basename>.disn_<isub>`.
  - **Format**: 
    - 2D ASCII: `I8, ',' , E20.12, ',' , E20.12` per line → `node_id, Ux, Uy` (comma-separated).
    - 3D ASCII: `I8, 3E20.12` → space-separated `node_id Ux Uy Uz`.
    - Binary: unformatted stream of `(Ux, Uy)` (2D) or `(Ux, Uy, Uz)` (3D) per node.
  - **Unit scaling**: Yes — divided by `1000.0` to mm when `*Key_Unit_System = 2`.

  Example (`Test_Examples/2d_fem/fem_2d.disn_1`):
  ~~~bash
         1,  0.000000000000E+00,  0.000000000000E+00
         2,  0.309792344217E-03,  0.200022243102E-02
         3,  0.306508276247E-04, -0.422207813393E-04
         4,  0.661986220995E-04, -0.714963121533E-04
         5,  0.116443200093E-03, -0.937367347085E-04
         6,  0.181700145801E-03, -0.103665261266E-03
         7,  0.259776862462E-03, -0.974530658974E-04
         8,  0.347292854805E-03, -0.718664920552E-04
         9,  0.440084114000E-03, -0.246253768389E-04
        10,  0.533547239550E-03,  0.455736716968E-04
  ~~~

  3D example (`Test_Examples/3d_uniaxial_tension/3d_uniaxial_tension.disn_1`):
  ~~~bash
         1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
         2 -0.158000824552E-02  0.000000000000E+00  0.000000000000E+00
         3 -0.599530622114E-04 -0.198353082554E-05  0.000000000000E+00
         4 -0.119662125493E-03 -0.473063658899E-05  0.000000000000E+00
         5 -0.178898410344E-03 -0.864250808527E-05  0.000000000000E+00
         6 -0.237351373110E-03 -0.138088347514E-04  0.000000000000E+00
         7 -0.294744765137E-03 -0.201261481292E-04  0.000000000000E+00
         8 -0.350861002901E-03 -0.273768250648E-04  0.000000000000E+00
         9 -0.405543687855E-03 -0.352738546285E-04  0.000000000000E+00
        10 -0.458697601816E-03 -0.434876069018E-04  0.000000000000E+00
  ~~~

+ ***<font color=green>.dien_<i></font>** - Enriched-DOF displacement vector (XFEM).

  - **Pattern**: `<basename>.dien_<isub>`.
  - **Format**: Same as `.disn_` (2D comma-separated, 3D space-separated) but values are **zero** for non-enriched nodes.

+ ***<font color=green>.edei_<i></font>** - Enriched-DOF index mapping. Format: `node_id, enriched_index` per row.

+ ***<font color=green>.disg_<i></font>** - Gauss-point displacements.

  - **Format**: 2D = `(I8, 2E20.12)` per row → `gp_id, Ux, Uy`. 3D adds `Uz`.

+ ***<font color=green>.dipc_<i></font>** - Cylindrical-coordinate raw DOF vector. Enabled when `*Key_CoorSys = 2`.
+ ***<font color=green>.dinc_<i></font>** - Cylindrical-coordinate nodal displacement, same indexing scheme as `.disn_`.

+ ***<font color=green>.veln_<i></font>** - Per-node velocity (dynamic analysis only).

  - **Format**: 2D = `(I8, ',', E20.12, ',', E20.12)` → `node_id, vx, vy`. 3D = `(I8, 3E20.12)` → `node_id, vx, vy, vz`.
  - **Units**: m/s (SI) or mm/s (mm-ton-s).

+ ***<font color=green>.acln_<i></font>** - Per-node acceleration (dynamic analysis only).

  - **Format**: identical layout to `.veln_`. Units: m/s² or mm/s².

<a id="section-Stress-files"></a>

## Stress files (.strn_, .strg_, .stnc_, .sttn_, .elss_)

+ ***<font color=green>.strn_<i></font>** - Nodal stress (the most-used output for cloud plots).

  - **Pattern**: `<basename>.strn_<isub>`.
  - **Format**:
    - 2D ASCII: `(I8, 4E20.12)` → `node_id, σxx, σyy, σxy, σvm`.
    - 3D ASCII: `(6E20.12)` → `σxx, σyy, σzz, σxy, σyz, σxz` (no node ID; ordering is implicit, one row per node).
    - Binary: unformatted stream.
  - **Units**: Pa (SI), or MPa (multiplied by `1.0e6` when `*Key_Unit_System = 2`).

  Example (`Test_Examples/2d_fem/fem_2d.strn_1`):
  ~~~bash
         1 -0.156417611316E+09 -0.521392037720E+09  0.132479407963E+09  0.517120498173E+09
         2 -0.314927052617E+09 -0.131351166901E+09 -0.151570950392E+09  0.379456717788E+09
         3  0.132437244468E+08 -0.434498039671E+09  0.491542972342E+08  0.449407109240E+09
         4 -0.355598495381E+08 -0.374743859893E+09  0.100279409204E+09  0.398170199644E+09
         5 -0.604493714661E+08 -0.352794869810E+09  0.137762121211E+09  0.404632332692E+09
         6 -0.100352136399E+09 -0.312552949135E+09  0.169284080617E+09  0.402946449247E+09
         7 -0.144585411303E+09 -0.264103676073E+09  0.187784543492E+09  0.397818122474E+09
         8 -0.189543490576E+09 -0.211017942701E+09  0.192461730240E+09  0.389336335768E+09
  ~~~

+ ***<font color=green>.strg_<i></font>** - Gauss-point stress.

  - **Format**:
    - 2D ASCII: `(I8, 4E20.12)` → `gp_id, σxx, σyy, σxy, σvm`.
    - 3D ASCII: `(I8, 7E20.12)` → `gp_id, σxx, σyy, σzz, σxy, σyz, σxz, σvm`.
    - Binary: stream of `(σxx, σyy, σxy, σvm)` per Gauss point (2D only).
  - **Units**: Pa or MPa.

+ ***<font color=green>.stnc_<i></font>** - Nodal stress in cylindrical coordinates (3D only). Format: `node_id, σrr, σθθ, σzz, σrθ, σθz, σrz, σvm`.

+ ***<font color=green>.sttn_<i></font>** - Thermal-stress component only at each node. Enabled by `*Key_Thermal_Stress = 1`. Same column layout as `.strn_`.

+ ***<font color=green>.elss_<i></font>** - Element stress state flag (1-3 stress condition).

  - **Format**: one `I2` integer per element. `0` = element does not satisfy the 1-3 condition; `1` = does.

<a id="section-Strain-files"></a>

## Strain files (.sran_, .srac_)

+ ***<font color=green>.sran_<i></font>** - Nodal strain tensor.

  - **Format**: 2D = `(I8, 4E20.12)` → `node_id, εxx, εyy, εxy, εvm`. 3D = `(I8, 7E20.12)` (adds εzz, εyz, εxz).
  - **Units**: dimensionless.

+ ***<font color=green>.srac_<i></font>** - Nodal strain in cylindrical coordinates (3D). Columns: `node_id, εrr, εθθ, εzz, εrθ, εθz, εrz, εvm`.

<a id="section-Gauss-point-files"></a>

## Gauss-point files (.gcor_, .disg_, .damg_, .elgn_, .strg_)

+ ***<font color=green>.gcor_<i></font>** - Gauss-point coordinates.

  - **Format**: 2D = `(E20.12)` per Gauss point alternating X, Y, X, Y, …. Equivalently: `(I8, 2E20.12)` → `gp_id, X, Y`.
  - **Units**: m (SI) or mm (mm-ton-s).

+ ***<font color=green>.disg_<i></font>** - Gauss-point displacements. See the `.disg_` entry above; same format.

+ ***<font color=green>.damg_<i></font>** - Per-Gauss-point damage factor.

  - **Format**: `(E20.12)` per Gauss point, one value per line.

+ ***<font color=green>.elgn_<i></font>** - Number of Gauss points per element.

  - **Format**: one `I10` per element, one value per line.

+ ***<font color=green>.strg_<i></font>** - See "Stress files" above.

<a id="section-SIF-fracture-energy-files"></a>

## SIF / fracture-energy files (.sifs_, .sift_, .prst_, .enee)

+ ***<font color=green>.sifs_<i></font>** - Stress intensity factors (KI, KII) for both tips of every crack.

  - **Pattern**: `<basename>.sifs_<isub>`.
  - **Format**: One row per crack, `(4E20.12)` → `KI_tip1, KII_tip1, KI_tip2, KII_tip2`.
  - **Units**: Pa·√m (SI) or MPa·√m (mm-ton-s, divided by `1.0e-6 * sqrt(1000)`).

  Example (`Test_Examples/2d_L_shaped_specimen/2d_L_shaped_specimen.sifs_1` — single crack):
  ~~~bash
    0.000000000000E+00  0.000000000000E+00  0.100294273438E+07  0.438159774968E+05
  ~~~

  Here crack 1 has KI₁ = 0, KII₁ = 0 at the first tip (a newly-formed tip), and KI₂ = 1.003e6 Pa·√m, KII₂ = 4.38e4 Pa·√m at the second (existing) tip.

+ ***<font color=green>.sift_<i></font>** - SIF time-history (NEWFTU-2026041702). Same format as `.sifs_` but accumulated per-step for plotting.
+ ***<font color=green>.prst_<i></font>** - Propagation speed and SIFs at a specified tip (NEWFTU-2026060901). Format: per-tip values, see `Save_SIFs_KI_and_KII.f90`.
+ ***<font color=green>.enee_<i></font>** - Fracture-energy increment per substep (used in cohesive-crack modelling).

<a id="section-Crack-coordinates"></a>

## Crack coordinates (.crax/.cray/.craz, .crxo/.cryo, .cnox/.cnoy/.cnoz)

+ ***<font color=green>.crax_<i></font>** - 2D / 3D X-coordinates of crack polyline.

  - **Format**: One row per crack, `(2000E20.12)` → all X-coordinates of that crack's polyline on a single line. For 3D, the row count equals `num_Crack`; each row contains `Each_Cr_Poi_Num(i)` reals.

+ ***<font color=green>.cray_<i></font>** - Y-coordinates, same format as `.crax_`.
+ ***<font color=green>.craz_<i></font>** - Z-coordinates (3D only). Same format.

  Example (`Test_Examples/2d_L_shaped_specimen/2d_L_shaped_specimen.crax_1`):
  ~~~bash
    0.137028414927E-04  0.247500000000E-01
  ~~~

  This row is the entire crack 1 polyline X-coordinates (2 points in this early step).

+ ***<font color=green>.crxo_<i>, .cryo_<i></font>** - Original (un-edge-disposed) crack coordinates, written when the solver is restarted from a previous step.

+ ***<font color=green>.cnox_<i>, .cnoy_<i>, .cnoz_<i></font>** - 3D crack-meshed-node coordinates. Format: one row per crack, `(50000E20.12)` per coordinate direction. Used for visualisation of the crack surface in 3D.

<a id="section-Crack-meshes"></a>

## Crack meshes (.cms1/.cms2/.cms3, .cmso, .cmap, .cape, .capf, .ctap)

+ ***<font color=green>.cms1_<i>, .cms2_<i>, .cms3_<i></font>** - Crack-mesh connectivity (3D crack surface triangulation).

  - **Format**: `I10` integers, one per line. Index into the crack-node coordinate arrays (`.cnox_`, `.cnoy_`, `.cnoz_`). Triangles connect `(cms1, cms2, cms3)` node triples.
  - **Indexing**: 1-based integers.

+ ***<font color=green>.cmso_<i></font>** - Crack-mesh outline (the boundary polygon of the crack surface).

+ ***<font color=green>.cmap_<i></font>** - Crack aperture at each meshed crack node (3D). Format: `(50000E20.12)` per crack.

+ ***<font color=green>.cape_<i></font>** - Crack aperture at each calculation point (2D and 3D). Format: `(2000E20.12)` per crack.

  - **Units**: m (SI) or mm (mm-ton-s).

+ ***<font color=green>.capf_<i></font>** - Aperture at fluid-element calculation points (3D HF). Same per-crack layout as `.cape_`.

+ ***<font color=green>.ctap_<i></font>** - Crack tangential (shear) opening, used in cohesive and contact analyses.

<a id="section-Crack-vertex-tangent-files"></a>

## Crack vertex and tangent files (.cvpx/.cvpy/.cvpz, .cvxx/.cvxy/.cvxz, ...)

These nine files describe the per-vertex coordinate frame of 3D cracks. Each is a per-crack row of `(50000E20.12)` reals, one file per coordinate direction.

+ ***<font color=green>.cvpx_<i>, .cvpy_<i>, .cvpz_<i></font>** - Vertex coordinates of the 3D crack outline.
+ ***<font color=green>.cvxx_<i>, .cvxy_<i>, .cvxz_<i></font>** - Local x-axis tangent vector at each vertex.
+ ***<font color=green>.cvyx_<i>, .cvyy_<i>, .cvyz_<i></font>** - Local y-axis tangent vector at each vertex.
+ ***<font color=green>.cvzx_<i>, .cvzy_<i>, .cvzz_<i></font>** - Local z-axis tangent vector at each vertex.

<a id="section-Crack-calculation-point-files"></a>

## Crack calculation-point files (.apex/.apey, .cori, .cpre, .cvel, .cqua, .cohx/.cohy, .crrd)

+ ***<font color=green>.apex_<i>, .apey_<i></font>** - Calculation-point coordinates (2D HF).

  - **Format**: per-crack row of `(2000E20.12)` reals. `.apex_` = X, `.apey_` = Y.
  - **Units**: m (SI) or mm (mm-ton-s).

+ ***<font color=green>.cori_<i></font>** - Calculation-point orientation angles (radians).

  - **Format**: per-crack row of `(2000E20.12)` reals.

+ ***<font color=green>.cpre_<i></font>** - Pore / fluid pressure at calculation points (HF only).

  - **Format**: per-crack row of `(2000E20.12)` reals. Units: Pa (SI), MPa (mm-ton-s).

+ ***<font color=green>.cvel_<i></font>** - Flow velocity at calculation points (HF). Units: m/s or mm/s.
+ ***<font color=green>.cqua_<i></font>** - Flow quantity (volumetric flow rate) at calculation points (HF). Units: m³/s.
+ ***<font color=green>.cohx_<i>, .cohy_<i></font>** - Cohesive tractions (only for cohesive-crack analyses).
+ ***<font color=green>.crrd_<i></font>** - Crack radius for radial cracks. Enabled by `*Key_Save_Crack_Radius = 1`.

<a id="section-XFEM-enrichment-files"></a>

## XFEM enrichment files (.ennd, .elty, .posi, .elts, .posh, .ennh, ...)

These describe the per-node / per-element XFEM enrichment topology. All are step-suffixed.

+ ***<font color=green>.ennd_<i></font>** - Enriched-node-type per crack.

  - **Format**: `Num_Node` rows × `num_Crack` columns in `(200I10)` format.
  - **Type codes (2D)**: `0` = no enrichment, `1` = tip-enriched, `2` = Heaviside, `3` = junction, `4` = Heaviside-tip, `6` = hole-junction.

+ ***<font color=green>.enns_<i></font>** - Enriched-node-type for cross-interface analyses.
+ ***<font color=green>.ennh_<i></font>** - Enriched-node-type for hole analyses. Format: `(200I10)`, columns = `num_Hole`.
+ ***<font color=green>.ennj_<i></font>** - Enriched-node-type for inclusion analyses. Columns = `num_Inclusion`.

+ ***<font color=green>.elty_<i></font>** - Element type per crack.

  - **Format**: `Num_Elem` rows × `num_Crack` columns in `(200I10)` format.
  - **Type codes (2D)**: `0` = none; `1` = tip; `2` = fully cracked (no kink); `3` = fully cracked (with kink); `4` = junction; `5` = complex junction; `6` = crack-hole junction.
  - **Type codes (3D)**: `0` = none; `1` = tip; `2` = Heaviside.

+ ***<font color=green>.elts_<i></font>** - Element type for cross-interface.
+ ***<font color=green>.elth_<i></font>** - Element type for holes. Format: `(200I10)`, columns = `num_Hole`.
+ ***<font color=green>.eltj_<i></font>** - Element type for inclusions.

  Example (`Test_Examples/2d_L_shaped_specimen/2d_L_shaped_specimen.elty_1` — first 10 rows):
  ~~~bash
           2
           0
           0
           0
           0
           0
           0
           0
           0
           0
  ~~~

  Each value is the integer type code for that element w.r.t. crack 1 (column-wise, since this is a single-crack problem).

+ ***<font color=green>.posi_<i>, .poss_<i>, .posh_<i>, .posj_<i></font>** - Signed-distance / position signs (per node per feature). Same `(200I10)` row layout as `.ennd_`. Values are `±1` (side of the discontinuity the node lies on).

+ ***<font color=green>.njel_<i></font>** - Crack-node to junction-element association. Integer pairs.

+ ***<font color=green>.nods_<i></font>** - Node-to-host-element mapping for cross-interface problems.

<a id="section-Crack-tip-baseline-files"></a>

## Crack-tip and baseline files (.celt, .celv, .celj, .celc, .ctty, .blab, .blvx/.blvy/.blvz, .tere)

These files describe the local geometry of enriched elements around the crack tip.

+ ***<font color=green>.celt_<i></font>** - Tip-element coordinates per element.
+ ***<font color=green>.celv_<i></font>** - Vertex coordinates per element.
+ ***<font color=green>.celj_<i></font>** - Junction coordinates per element per crack.
+ ***<font color=green>.celc_<i></font>** - Element crack-intersection coordinates.
+ ***<font color=green>.ctty_<i></font>** - Crack-tip type per crack.
+ ***<font color=green>.blab_<i></font>** - Crack baseline vectors. Format: 6 values per element (two endpoint coordinates X, Y, Z).
+ ***<font color=green>.blvx_<i>, .blvy_<i>, .blvz_<i></font>** - Baseline direction-vector components (one file per axis).
+ ***<font color=green>.tere_<i></font>** - Tip-enriched node element number. Format: `enriched_node_number, crack_id`.

<a id="section-Crack-node-S1-vector-files"></a>

## Crack-node S1 vector files (.cndx/.cndy/.cndz)

+ ***<font color=green>.cndx_<i>, .cndy_<i>, .cndz_<i></font>** - Crack-node S1 vector components used in the CFCP-2 crack-propagation criterion (maximum principal tensile stress). Enabled when `*CFCP = 2`.

<a id="section-Hole-inclusion-files"></a>

## Hole and inclusion files (.hlcr, .ehcr, .jzcr, .jzpx/.jzpy)

Holes and inclusions have **no step suffix** — written once per analysis.

+ ***<font color=green>.hlcr</font>** - Circular hole coordinates.

  - **Format**: `num_Circ_Hole` rows, `(3E20.12)` per row → `x_center, y_center, z_center`.

+ ***<font color=green>.ehcr</font>** - Elliptical hole parameters.

  - **Format**: `num_Ellip_Hole` rows, `(5E20.12)` per row → `x_c, y_c, a, b, θ` (center, semi-axes, rotation angle).

+ ***<font color=green>.jzcr</font>** - Circular inclusion coordinates.

  - **Format**: `num_Circ_Incls` rows, `(3E20.12)` per row → `x, y, r` (center and radius).

+ ***<font color=green>.jzpx, .jzpy</font>** - Polygonal inclusion vertex coordinates (X and Y). Format: per-inclusion rows of `(2000E20.12)`.

<a id="section-Cross-interface-files-cscr"></a>

## Cross-interface files (.cscr)

+ ***<font color=green>.cscr</font>** - Cross-interface point coordinates.

  - **Format**: per-cross row of `(2E20.12)` → `x, y`.

<a id="section-Natural-fracture-files"></a>

## Natural-fracture files (.nfcx/.nfcy/.nfcz, .ncrx/.ncry)

+ ***<font color=green>.nfcx, .nfcy, .nfcz</font>** - Natural-fracture vertex coordinates (3D). Format: per-fracture row of `(50000E20.12)` per axis.

+ ***<font color=green>.ncrx, .ncry</font>** - Natural-crack endpoint coordinates (2D). Format: `x1, x2` comma-separated per crack.

<a id="section-Energy-history-files"></a>

## Energy history files (.ener, .edye, .edtm, .edap)

+ ***<font color=green>.ener</font>** - Energy balance per load step. **Append-mode** — written once per step.

  - **Header** (written at `iter = 1`):
    ~~~bash
    # iter        Elastic_Strain_Energy Fracture_Energy       External_Work         Residual_Energy      Normalized_Residual    Crack_Propagation_Length
    ~~~
  - **Data rows** (`(I10, 6(2X, E20.12))`):
    `iter, Elastic_Strain_Energy, Fracture_Energy, External_Work, Residual_Energy, Normalized_Residual_Energy, Crack_Propagation_Length`
  - **Units**: J (SI), or mJ (mm-ton-s); length in m or mm.

  Example (`Test_Examples/2d_L_shaped_specimen/2d_L_shaped_specimen.ener`):
  ~~~bash
  # iter        Elastic_Strain_Energy Fracture_Energy       External_Work         Residual_Energy      Normalized_Residual    Crack_Propagation_Length
           1    0.472931017820E+01   -0.129866029815E-15    0.473748321243E+01    0.817303423107E-02    0.172518484279E+00   -0.346944695195E-17
           2    0.403905163275E+01    0.935783374928E+00    0.498766375932E+01    0.128287516428E-01    0.257209632844E+00    0.250000002500E-01
           3    0.349579005177E+01    0.187156673676E+01    0.541280434455E+01    0.454475560179E-01    0.839630496966E+00    0.500000001500E-01
           4    0.299507968076E+01    0.280735009297E+01    0.588459158372E+01    0.821618099953E-01    0.139621941177E+01    0.749999999000E-01
           5    0.252293261486E+01    0.374313345230E+01    0.635674090155E+01    0.906748343931E-01    0.142643590163E+01    0.999999997333E-01
           6    0.200212858936E+01    0.467891680596E+01    0.682556546132E+01    0.144520066004E+00    0.211733470029E+01    0.124999999415E+00
           7    0.155959488881E+01    0.561470017729E+01    0.732121007053E+01    0.146915004425E+00    0.200670385100E+01    0.149999999569E+00
           8    0.114477972662E+01    0.655048354286E+01    0.788939019863E+01    0.194126929148E+00    0.246060752809E+01    0.174999999569E+00
  ~~~

  Reading: at step 8, the elastic-strain energy has dropped to 1.14 J while fracture-energy has accumulated to 6.55 J, with crack propagating a total of 0.175 m.

+ ***<font color=green>.edye</font>** - Explicit-dynamics energy history.

  - **Header**: `# iter           time               kinetic_energy (T)     internal_energy (U)    total_energy (T+U)   external_Work (W)`.
  - **Data**: `iter, time, T, U, T+U, W` (`(I10, 5E20.12)`).

+ ***<font color=green>.edtm</font>** - Explicit-dynamics time history. Format: `(I10, E20.12)` → `iter, time`.

+ ***<font color=green>.edap</font>** - Explicit-dynamics max/min/avg aperture per crack (3D). Format: `(E20.12, I10, 3E20.12)` → `time, crack_id, max_aperture, min_aperture, avg_aperture`.

<a id="section-Crack-length-history"></a>

## Crack length history (.dcrl, .dprl, .idtm, .iite, .ilth, .ipre)

+ ***<font color=green>.dcrl</font>** - Total crack length history (dynamic 2D).

  - **Header**: `iter        c_Time(s)            Crack_1_Length(m)`.
  - **Data**: `(I10, 2E20.12)` → `iter, time, total_crack_length`.

+ ***<font color=green>.dprl</font>** - Crack propagation length history (the incremental growth, **excluding** the initial crack length).

  - **Format**: same as `.dcrl`. Use this for the "growth length" vs energy curves.

<a id="section-Dynamic-files"></a>

## Dynamic files (.idtm, .edtm, .dcrl, .dprl)

Implicit-dynamics and explicit-dynamics share several time-history files:

+ ***<font color=green>.idtm</font>** - Implicit-dynamics time history. Format: `iter, time` (`(I10, E20.12)`).

The other dynamic files (`.edtm`, `.dcrl`, `.dprl`) are covered in their respective sections.

<a id="section-Lumped-mass-files"></a>

## Lumped mass files (.lpmf_, .lpmx_, .lpma_)

+ ***<font color=green>.lpmf_<i></font>** - Lumped mass at each node from the standard FEM mass matrix.

  - **Format**: one `(E20.12)` per node.

+ ***<font color=green>.lpmx_<i></font>** - Lumped mass from the enriched (XFEM) mass matrix. Zero for non-enriched nodes.

+ ***<font color=green>.lpma_<i></font>** - Superposition (FEM + enriched) lumped mass at each node.

<a id="section-HF-output-files"></a>

## Hydraulic-fracturing output files (.hftm, .injp, .wbfp, .wbpt, .ihft, .icpt, .icpt_LS, .icpt_CT)

+ ***<font color=green>.hftm</font>** - HF analysis time-history log. **Append-mode**.

  - **Header** (at first call `(1, 1, 1)`):
    ~~~bash
        imf   |   ifra   | total_ter|   time   
    ~~~
  - **Data**: `(3I10, 1F18.5)` → `imf, ifra, Counter_Iter, c_Time`.
  - **Cleared** at the start of each new analysis run by deleting the file (lines 36-40 of `Save_HF_time.f90`).

+ ***<font color=green>.injp</font>** - HF injection-pressure history. **Append-mode**.

  - **Format**: `(2E20.12)` per row → `time, pressure`.
  - **Units**: seconds and Pa (SI), or seconds and MPa (mm-ton-s).

+ ***<font color=green>.wbfp</font>** - Wellbore fracturing path (XYZ coordinates of wellbore stages).

  - **Format**: 3 reals per wellbore point (XYZ), per stage.

+ ***<font color=green>.wbpt</font>** - Wellbore/stage/proppant pressure vs time.

  - **Format**: 5 columns → `i_WB, i_Stage, i_Prop, time, pressure`.

+ ***<font color=green>.ihft</font>** - HF time per fracture step.

  - **Format**: per-step `(F12.5)` row (only positive entries written).
  - **Header**: first row is the string `HF time of each time step`.

+ ***<font color=green>.icpt</font>** - CPU time per step (total wall time).
+ ***<font color=green>.icpt_LS</font>** - CPU time spent in the linear solver per step.
+ ***<font color=green>.icpt_CT</font>** - CPU time spent in cohesive/contact computations per step.

+ ***<font color=green>.iite</font>** - Number of iterations per fracture step.
+ ***<font color=green>.ilth</font>** - Crack-1 length per fracture step.
+ ***<font color=green>.ipre</font>** - Maximum pressure in crack 1 per fracture step.

<a id="section-Proppant-output-files"></a>

## Proppant output files (.ccon_, .pokf_, .cond_, .wpnp_, .epcr_, .Saved_Filename)

All proppant files are enabled by `*Key_Propp_Trans = 1`.

+ ***<font color=green>.ccon_<i></font>** - Proppant concentration at calculation points.

+ ***<font color=green>.pokf_<i></font>** - Proppant permeability (`kf = Cf / w`, concentration divided by aperture).

+ ***<font color=green>.cond_<i></font>** - Proppant conductivity.

+ ***<font color=green>.wpnp_<i></font>** - Propped aperture at zero closure.

+ ***<font color=green>.epcr_<i></font>** - Element proppant coordinates.

  - **Format**: per-element row `(I8, 4E20.12)` → `element_id, x, y, z, w` (centroid X, Y, Z and proppant mass `w`).
  - **Only** elements with `Elem_Proppant_Coor(i,1) = 1` are written.

+ ***<font color=green>.Saved_Filename</font>** - Flow rate at a specific wellbore point. The filename is taken from a keyword-supplied name (not auto-generated).

  - **Format**: `ifra, time, flow_Q` per row.

<a id="section-Contact-cohesive-element-state"></a>

## Contact, cohesive and element-state files (.elcs_, .elco_, .elss_, .kiel_)

+ ***<font color=green>.elcs_<i></font>** - Contact state per element. Integer per element:

  - `0` = open
  - `1` = stick
  - `2` = slip

+ ***<font color=green>.elco_<i></font>** - Cohesive state per element. Integer per element: `0` = inactive, `1` = active.

+ ***<font color=green>.elss_<i></font>** - Element 1-3 stress-condition flag. Integer per element: `0` = not satisfied, `1` = satisfied. See `.elss_` in Stress files above.

+ ***<font color=green>.kiel_<i></font>** - Killed (broken / deactivated) element IDs.

  - **Format**: one `I10` integer per line; accumulated from all prior steps plus current step.

<a id="section-DOF-surface-load-files"></a>

## DOF / surface-load force files (.fxdf_, .fydf_, .fzdf_, .fxsl_, .fysl_, .fzsl_)

+ ***<font color=green>.fxdf_<i></font>** - X-DOF internal force vector.

  - **Format**: one `(E20.12)` per node.

+ ***<font color=green>.fydf_<i></font>** - Y-DOF internal force.
+ ***<font color=green>.fzdf_<i></font>** - Z-DOF internal force (3D only).

+ ***<font color=green>.fxsl_<iSL>_<i></font>** - Surface-load X-component per surface load `iSL`. Same format.
+ ***<font color=green>.fysl_<iSL>_<i></font>** - Surface-load Y-component per `iSL`.
+ ***<font color=green>.fzsl_<iSL>_<i></font>** - Surface-load Z-component per `iSL` (3D only).

<a id="section-Stiffness-matrix-files"></a>

## Stiffness-matrix files (.skxf, .csrn_, .csra_, .csrj_, .csri_)

+ ***<font color=green>.skxf</font>** - XFEM element stiffness matrices (binary). Each XFEM element's flattened K-matrix is written as a binary record. Async-write mode is used for performance.

+ ***<font color=green>.csrn_<i></font>** - CSR matrix dimensions.

  - **Format**: 3 lines, each `(I0)`:
    ~~~bash
    K_CSR_NNZ
    num_FreeD
    XFEM_Start_DOF
    ~~~

+ ***<font color=green>.csra_<i></font>** - CSR value array (`K_CSR_aa`). Format: one `(E25.16E3)` value per line, `K_CSR_NNZ` rows.

+ ***<font color=green>.csrj_<i></font>** - CSR column-index array (`K_CSR_ja`). Format: one `(I0)` value per line, `K_CSR_NNZ` rows.

+ ***<font color=green>.csri_<i></font>** - CSR row-pointer array (`K_CSR_ia`). Format: one `(I0)` value per line, `num_FreeD + 1` rows.

  These four CSR files together define the assembled CSR sparse stiffness matrix. They can be re-loaded by external linear-system codes for further processing.

<a id="section-3D-fluid-crack-files"></a>

## 3D fluid-crack files (.fenx/.feny/.fenz, .fnnx/.fnny/.fnnz, .fnux/.fnuy/.fnuz, .fnlx/.fnly/.fnlz, .fexx/.fexy/.fexz/.feyx/.feyy/.feyz/.fezx/.fezy/.fezz, .cpfn_, .cpno_, .ccpx/.ccpy/.ccpz, .cnlx/.cnly/.cnlz, .cmse_)

3D hydraulic-fracturing uses fluid elements embedded in the crack surface. All files here are step-suffixed and use `(E20.12)` per row or `(50000E20.12)` per crack.

+ ***<font color=green>.fenx_<i>, .feny_<i>, .fenz_<i></font>** - Fluid-element normal vectors (per-element).
+ ***<font color=green>.fnnx_<i>, .fnny_<i>, .fnnz_<i></font>** - Fluid-node normal vectors (per-node).
+ ***<font color=green>.fnux_<i>, .fnuy_<i>, .fnuz_<i></font>** - Up-face displacement vectors at crack nodes.
+ ***<font color=green>.fnlx_<i>, .fnly_<i>, .fnlz_<i></font>** - Low-face displacement vectors at crack nodes.
+ ***<font color=green>.fexx_<i>, .fexy_<i>, .fexz_<i>, .feyx_<i>, .feyy_<i>, .feyz_<i>, .fezx_<i>, .fezy_<i>, .fezz_<i></font>** - Fluid-element local-coordinate axes (x, y, z basis at each element centroid).
+ ***<font color=green>.cpfn_<i></font>** - Fluid node count per crack (one integer per crack).
+ ***<font color=green>.cpno_<i></font>** - Fluid node numbers per crack.
+ ***<font color=green>.ccpx_<i>, .ccpy_<i>, .ccpz_<i></font>** - Fluid-node coordinates (X, Y, Z).
+ ***<font color=green>.cnlx_<i>, .cnly_<i>, .cnlz_<i></font>** - Crack-node local coordinates (in the fluid-element frame).
+ ***<font color=green>.cmse_<i></font>** - Element containing each crack node.

<a id="section-Miscellaneous-files"></a>

## Miscellaneous files (.fraz, .seed, .post, .fdcu, .fccu, .rbco_, .fdvl, .ecfv, .sccx/.sccy/.sccz, .scdx/.scdy/.scdz)

+ ***<font color=green>.fraz</font>** - Fracture-zone bounding box.

  - **Format 2D**: `F12.5` × 4 → `min_X, max_X, min_Y, max_Y`.
  - **Format 3D**: `F12.5` × 6 → `min_X, max_X, min_Y, max_Y, min_Z, max_Z`.

+ ***<font color=green>.seed</font>** - Random-number seed (used by stochastic analyses). Format: plain text `I11`.

+ ***<font color=green>.post</font>** - Header for the legacy MATLAB post-processor.

  - **Format**: single line with six values:
    `Key_Analysis_Type, Key_TipEnrich, Key_Data_Format, Key_Heaviside_Value, Key_Hole_Value, Ave_Elem_L`

+ ***<font color=green>.fdcu</font>** - Force-displacement curve at a specified node.

  - **Header**: `node |  isub | lambda | dis_of_node`.
  - **Data**: `node_id, isub, λ, displacement`.

+ ***<font color=green>.fccu</font>** - Force-COD (crack-opening-displacement) curve.

  - **Header**: `crack  |  isub  |   lambda   |   COD`.

+ ***<font color=green>.rbco_<i></font>** - Rigid-ball (circle) coordinates (for rigid-ball contact problems). Format: `x, y, radius` per ball.

+ ***<font color=green>.fdvl_<i></font>** - Fluid-velocity per-node (for HF post-processing). Format: 2D = `node_id, vx, vy`; 3D = `node_id, vx, vy, vz`.

+ ***<font color=green>.ecfv_<i></font>** - Element-centroid field value (used by some HF post-processors).

+ ***<font color=green>.sccx_<i>, .sccy_<i></font>** - 2D stress-corrosion-cracking X / Y coordinates.
+ ***<font color=green>.scdx_<i>, .scdy_<i></font>** - 2D stress-corrosion displacement X / Y.

<a id="section-Console-log-files"></a>

## Console / log files (PhiPsi_Console_Window.log, current_folder.dat)

+ ***<font color=green>PhiPsi_Console_Window.log</font>** - Console-window text capture. Written next to the executable when `*Key_Save_Window_Log = 1` (or equivalent).

+ ***<font color=green>current_folder.dat</font>** - Single-line text file containing the last `Work_Directory`. Used by PPView to remember the most-recent folder. Written at `<PhiPsi current directory>\current_folder.dat` and `<PhiPsi current directory>\python_tools\current_folder.dat`.

---

# PPView file-reading map

PPView (`PPView/PPView/`) consumes the output files above and renders them as cloud plots, curves, vectors, and 3D meshes. This section lists the file-to-feature mapping discovered in `pp_utilities.py`, `pp_plotting.py`, and `child_window.py`.

<a id="section-Cloud-plot-rendering"></a>

## Cloud plot / contour rendering

PPView's cloud-plot window reads displacement / stress / strain fields and renders them as a colour-mapped surface over the mesh.

| Feature | Files read | Columns / fields used |
|---|---|---|
| Displacement magnitude | `*.disn_<i>` | `node_id, Ux, Uy` (2D comma) or `Ux, Uy, Uz` (3D space). Magnitude = `sqrt(Ux² + Uy² + …)`. |
| Displacement X / Y / Z | `*.disn_<i>` | single component. |
| Stress σxx / σyy / σxy / σvm | `*.strn_<i>` | 2D = 4 stress columns; 3D = 6 stress columns. |
| Stress σzz / σyz / σxz | `*.strn_<i>` | 3D-only columns. |
| Cylindrical stress σrr / σθθ / σrθ | `*.stnc_<i>` | 3D-only columns. |
| Strain | `*.sran_<i>` | 2D = 4 strain columns; 3D = 7. |
| Stress / strain at Gauss points | `*.strg_<i>`, `*.disg_<i>` | one row per Gauss point. |
| Enriched-node indicator | `*.ennd_<i>`, `*.elty_<i>` | integer per node / element. |
| Enriched-DOF displacement | `*.dien_<i>` | zero for non-enriched. |
| Full cloud (any field) | `*.vtk` | `SCALARS` fields written by `Save_vtk_file.f90` (stress, displacement, material ID, etc.). |
| Crack displacement (3D) | `*.fnux_<i>`, `*.fnlx_<i>` | up- and low-face displacement at crack nodes. |

<a id="section-Curve-plot-rendering"></a>

## Curve plot / line chart rendering

The curve-plot window (e.g. `Plot_Curve_Settings_Window` in `child_window.py`) reads time-history files and plots one column versus another.

| Feature | Files read | X-axis options | Y-axis options |
|---|---|---|---|
| Energy vs step / growth length | `*.ener` | `iter` (step #), or `Crack_Propagation_Length` (col 7) | `Elastic_Strain_Energy`, `Fracture_Energy`, `External_Work`, `Residual_Energy`, `Normalized_Residual_Energy` |
| Explicit-dynamics energy | `*.edye` | `time`, `iter` | `T`, `U`, `T+U`, `W` |
| SIF history | `*.sifs_<i>` (NEWFTU-2026041702 stack) | step # | `KI_tip1`, `KII_tip1`, `KI_tip2`, `KII_tip2` |
| Crack length / growth | `*.dcrl`, `*.dprl` | `iter`, `time` | `Crack_1_Length` |
| Displacement at a node | `*.disn_<i>` | node index | `Ux`, `Uy` |
| Velocity / acceleration at a node | `*.veln_<i>`, `*.acln_<i>` (NEWFTU-2026041702) | node index | single component |
| Lumped-mass | `*.lpmf_<i>`, `*.lpma_<i>` | node index | `mass` |

<a id="section-Vector-plot-rendering"></a>

## Vector plot / arrow rendering

The vector-plot window reads vector-valued fields and renders them as 3D arrows (via PyVista).

| Feature | Files read | Columns / fields used |
|---|---|---|
| Displacement arrows | `*.disn_<i>` | `(Ux, Uy)` (2D) or `(Ux, Uy, Uz)` (3D). |
| Velocity arrows | `*.veln_<i>` | as above. |
| Block / baseline vectors | `*.blvx_<i>`, `*.blvy_<i>`, `*.blvz_<i>` | one axis per file. |
| Fluid-element local axes | `*.fexx_<i>` … `*.fezz_<i>` | nine-component basis. |
| Original boundary forces | `*.focx`, `*.focy`, `*.focz` (input) | `node_id, force_value`. |

<a id="section-Crack-visualization"></a>

## Crack visualization (2D and 3D)

PPView's crack-overlay window reads crack geometry files and overlays them on the deformed mesh.

| Feature | Files read | Notes |
|---|---|---|
| 2D crack polyline | `*.crax_<i>`, `*.cray_<i>` | per-crack rows of X, Y coordinates. |
| 2D crack aperture | `*.cape_<i>` | aperture colour map. |
| 2D crack orientation | `*.cori_<i>` | used for tip-extension direction. |
| 2D enriched element types | `*.elty_<i>`, `*.ennd_<i>` | colour-code tip / Heaviside / junction elements. |
| 2D calc-point overlay | `*.apex_<i>`, `*.apey_<i>` | HF analysis pressure / aperture. |
| 3D crack surface mesh | `*.cnox_<i>`, `*.cnoy_<i>`, `*.cnoz_<i>` + `*.cms1_<i>`, `*.cms2_<i>`, `*.cms3_<i>` | triangle vertex coords + triangle connectivity. |
| 3D crack aperture | `*.cmap_<i>`, `*.cape_<i>` | colour map on the crack surface. |
| 3D crack outline | `*.cmso_<i>` | boundary polygon. |
| 3D natural fractures | `*.nfcx`, `*.nfcy`, `*.nfcz` | read once. |
| Clipped crack variants | `*.cnxc_<i>`, `*.cnyc_<i>`, `*.cnzc_<i>`, `*.cmc1_<i>` … | when the model is sliced by a plane. |
| Crack-aperture history | `*.cape_<i>` for `i = 1, 2, …, n` | frame-by-frame aperture animation. |

<a id="section-Animation-playback"></a>

## Animation playback

PPView generates GIF / PNG animations by reading sequential step-suffixed files.

| Source | Behavior |
|---|---|
| `<name>_NNNNN.vtk` sequence | Loaded into PyVista; each frame is a `vtkUnstructuredGrid` rendered in turn. |
| `<name>_<i>.disn_<i>` / `.strn_<i>` / `.crax_<i>` etc. sequence | Each frame's results are mapped onto the same `.node` / `.elem` mesh, then exported as PNG and stitched into a GIF (`WorkerThread_Generate_Animation` in `pp_utilities.py`). |
| 3D crack frame | `*.crax_<i>` / `*.cray_<i>` / `*.craz_<i>` for 2D; `*.cms1/2/3_<i>` for 3D. |

Animation parameters (`animation_interval`, `gif_scale_factor`) are configured in `child_window.py`'s animation tab.

---

# Data format examples

The next four sections show real excerpts from `Test_Examples/`. They are intentionally short (5–10 lines per file) so the manual stays readable; the corresponding full files in the repository are far larger.

<a id="section-Example-A-2d-fem"></a>

## Example A — 2D simple FEM (Test_Examples/2d_fem)

A quasi-static plane-stress FEM run on a square plate.

**`2d_fem.kpp`** — the keyword input file:
~~~bash
% Working directory.
*Work_Directory
X:\PhiPsi_Project\Test_Examples\2d_fem

% Filename of input files.
*Filename
fem_2d

% Quasi-static analysis.
*Key_Analysis_Type
1

% Plane stress.
*Key_Type_2D
1

% MUMPS solver.
*Key_SLOE
6

% One sub-step.
*Num_Substeps
1

% Material(1-E,2-v,3-density,4-thick,5-St,6-KIc,7-Sc,8-20(blank))
*Material_Para_1
70.0e9,0.3,2700.0,1.0,1.0e6,1.0e6,100.0e6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

% Save VTK files.
*Key_Save_vtk
1
~~~

**`fem_2d.node`** — 2D node coordinates (`X Y`):
~~~bash
  0.50000000E-01   0.00000000E+00
 -0.40000000E-01   0.30000000E-01
  0.49610478E-01   0.62290027E-02
  0.48447981E-01   0.12360952E-01
  0.46530622E-01   0.18300307E-01
  0.43888275E-01   0.23954527E-01
  0.40562109E-01   0.29235514E-01
  0.36603951E-01   0.34060986E-01
  0.32075470E-01   0.38355759E-01
  0.27047226E-01   0.42052914E-01
~~~

**`fem_2d.elem`** — 4-node quad elements (`n1 n2 n3 n4 mat`):
~~~bash
        1.       3.      81.      62.       1.
        3.       4.     100.      81.       1.
        4.       5.     119.     100.       1.
        5.       6.     138.     119.       1.
        6.       7.     157.     138.       1.
        7.       8.     176.     157.       1.
        8.       9.     195.     176.       1.
        9.      10.     214.     195.       1.
       10.      11.     233.     214.       1.
       11.      12.     252.     233.       1.
~~~

**`fem_2d.boux`** — left-edge nodes fixed in X:
~~~bash
        1.
       42.
       62.
       63.
       64.
       65.
       66.
       67.
       68.
       69.
~~~

**`fem_2d.focy`** — nodal Y-forces (downward traction on the top edge):
~~~bash
       22.   0.22500000E+05
      462.   0.22500000E+05
      482.   0.45000000E+05
      483.   0.45000000E+05
      484.   0.45000000E+05
      485.   0.45000000E+05
      486.   0.45000000E+05
      487.   0.45000000E+05
      488.   0.45000000E+05
      489.   0.45000000E+05
~~~

**`fem_2d.disn_1`** — output: nodal displacements (`node_id, Ux, Uy`, comma-separated):
~~~bash
         1,  0.000000000000E+00,  0.000000000000E+00
         2,  0.309792344217E-03,  0.200022243102E-02
         3,  0.306508276247E-04, -0.422207813393E-04
         4,  0.661986220995E-04, -0.714963121533E-04
         5,  0.116443200093E-03, -0.937367347085E-04
         6,  0.181700145801E-03, -0.103665261266E-03
         7,  0.259776862462E-03, -0.974530658974E-04
         8,  0.347292854805E-03, -0.718664920552E-04
         9,  0.440084114000E-03, -0.246253768389E-04
        10,  0.533547239550E-03,  0.455736716968E-04
~~~

**`fem_2d.strn_1`** — output: nodal stresses (`node_id, σxx, σyy, σxy, σvm`):
~~~bash
         1 -0.156417611316E+09 -0.521392037720E+09  0.132479407963E+09  0.517120498173E+09
         2 -0.314927052617E+09 -0.131351166901E+09 -0.151570950392E+09  0.379456717788E+09
         3  0.132437244468E+08 -0.434498039671E+09  0.491542972342E+08  0.449407109240E+09
         4 -0.355598495381E+08 -0.374743859893E+09  0.100279409204E+09  0.398170199644E+09
         5 -0.604493714661E+08 -0.352794869810E+09  0.137762121211E+09  0.404632332692E+09
         6 -0.100352136399E+09 -0.312552949135E+09  0.169284080617E+09  0.402946449247E+09
         7 -0.144585411303E+09 -0.264103676073E+09  0.187784543492E+09  0.397818122474E+09
         8 -0.189543490576E+09 -0.211017942701E+09  0.192461730240E+09  0.389336335768E+09
~~~

**`fem_2d_00001.vtk`** — output: ParaView-compatible VTK snapshot (header + first POINTS):
~~~bash
# vtk DataFile Version 4.0
X:\PhiPsi_Project\Test_Examples\2d_fem\fem_2d  - results from increment 00001
ASCII
DATASET UNSTRUCTURED_GRID

POINTS 861 double
    0.050000    0.000000    0.000000
   -0.040000    0.030000    0.000000
    0.049610    0.006229    0.000000
    0.048448    0.012361    0.000000
    0.046531    0.018300    0.000000
    0.043888    0.023955    0.000000
    0.040562    0.029236    0.000000
    0.036604    0.034061    0.000000
    0.032075    0.038356    0.000000
    0.027047    0.042053    0.000000
~~~

<a id="section-Example-B-2d-L-shaped-specimen"></a>

## Example B — 2D crack propagation (Test_Examples/2d_L_shaped_specimen)

A quasi-static plane-stress L-shaped panel with an initial crack propagating over 8 sub-steps. This example exercises nearly every output format.

**`2d_L_shaped_specimen.kpp`** (key lines):
~~~bash
*Work_Directory
X:\PhiPsi_Project\Test_Examples\2d_L_shaped_specimen

*Filename
2d_L_shaped_specimen

*Key_Analysis_Type
1                % Quasi-static.

*Key_Type_2D
1                % Plane stress.

*Key_SLOE
6                % MUMPS solver.

*Num_Substeps
8                % 8 propagation steps.

*Num_Crack
1

*CRACK2D_Coor_1
0.137028414927E-04,0.247500000000E-01

*Material_Para_1
20.0e9,0.3,2000.0,1.0,0.1e6,1.0e6,100.0e6,...

*Key_Save_vtk
1
~~~

**`2d_L_shaped_specimen.elem`** — first 10 rows:
~~~bash
        1.      60.      61.       3.       1.
       60.      59.      75.      61.       1.
       59.      58.      89.      75.       1.
       58.      57.     103.      89.       1.
       57.      56.     117.     103.       1.
       56.      55.     131.     117.       1.
       55.      54.     145.     131.       1.
       54.      53.     159.     145.       1.
       53.      52.     173.     159.       1.
       52.      51.     187.     173.       1.
~~~

**`2d_L_shaped_specimen.focy`** — applied Y-forces:
~~~bash
    512.   0.20000000E+04
    741.   0.10000000E+04
    742.   0.20000000E+04
    743.   0.20000000E+04
    744.   0.20000000E+04
~~~

**`2d_L_shaped_specimen.sifs_1`** — SIFs at step 1 (one crack, two tips):
~~~bash
  0.000000000000E+00  0.000000000000E+00  0.100294273438E+07  0.438159774968E+05
~~~

**`2d_L_shaped_specimen.crax_1`** — crack-1 X-coordinates (polyline row):
~~~bash
  0.137028414927E-04  0.247500000000E-01
~~~

**`2d_L_shaped_specimen.elty_1`** — element type per crack (first 10 elements):
~~~bash
         2
         0
         0
         0
         0
         0
         0
         0
         0
         0
~~~

**`2d_L_shaped_specimen.ener`** — full 8-step energy balance:
~~~bash
# iter        Elastic_Strain_Energy Fracture_Energy       External_Work         Residual_Energy      Normalized_Residual    Crack_Propagation_Length
         1    0.472931017820E+01   -0.129866029815E-15    0.473748321243E+01    0.817303423107E-02    0.172518484279E+00   -0.346944695195E-17
         2    0.403905163275E+01    0.935783374928E+00    0.498766375932E+01    0.128287516428E-01    0.257209632844E+00    0.250000002500E-01
         3    0.349579005177E+01    0.187156673676E+01    0.541280434455E+01    0.454475560179E-01    0.839630496966E+00    0.500000001500E-01
         4    0.299507968076E+01    0.280735009297E+01    0.588459158372E+01    0.821618099953E-01    0.139621941177E+01    0.749999999000E-01
         5    0.252293261486E+01    0.374313345230E+01    0.635674090155E+01    0.906748343931E-01    0.142643590163E+01    0.999999997333E-01
         6    0.200212858936E+01    0.467891680596E+01    0.682556546132E+01    0.144520066004E+00    0.211733470029E+01    0.124999999415E+00
         7    0.155959488881E+01    0.561470017729E+01    0.732121007053E+01    0.146915004425E+00    0.200670385100E+01    0.149999999569E+00
         8    0.114477972662E+01    0.655048354286E+01    0.788939019863E+01    0.194126929148E+00    0.246060752809E+01    0.174999999569E+00
~~~

<a id="section-Example-C-3d-uniaxial-tension"></a>

## Example C — 3D uniaxial tension (Test_Examples/3d_uniaxial_tension)

A 3D XFEM run on a rectangular block with one penny-shaped crack.

**`3d_uniaxial_tension.kpp`** (key lines):
~~~bash
*Work_Directory
X:\PhiPsi_Project\Test_Examples\3d_uniaxial_tension

*Filename
3d_uniaxial_tension

*Key_Dimension
3                % 3D problem.

*Key_Analysis_Type
1                % Quasi-static.

*Num_Crack
1

*CRACK3D_CIR_COOR_1
5.5,5.5,7.5,0,0,1,2.5

*Material_Para_1
20.0e9,0.3,2000.0,1.0,0.1e6,1.0e6,100.0e6,...

*Key_Save_vtk
1
~~~

**`3d_uniaxial_tension.node`** — 3D node coordinates (`X Y Z`):
~~~bash
  0.00000000E+00   0.00000000E+00   0.00000000E+00
  0.11000000E+02   0.00000000E+00   0.00000000E+00
  0.36666667E+00   0.00000000E+00   0.00000000E+00
  0.73333333E+00   0.00000000E+00   0.00000000E+00
  0.11000000E+01   0.00000000E+00   0.00000000E+00
  0.14666667E+01   0.00000000E+00   0.00000000E+00
  0.18333333E+01   0.00000000E+00   0.00000000E+00
  0.22000000E+01   0.00000000E+00   0.00000000E+00
~~~

**`3d_uniaxial_tension.elem`** — 8-node hexahedra (`n1 n2 n3 n4 n5 n6 n7 n8 mat`):
~~~bash
        1.       3.     121.     120.    1072.    1073.    6723.    4722.       1.
        3.       4.     150.     121.    1073.    1113.    7883.    6723.       1.
        4.       5.     179.     150.    1113.    1153.    9043.    7883.       1.
        5.       6.     208.     179.    1153.    1193.   10203.    9043.       1.
        6.       7.     237.     208.    1193.    1233.   11363.   10203.       1.
        7.       8.     266.     237.    1233.    1273.   12523.   11363.       1.
        8.       9.     295.     266.    1273.    1313.   13683.   12523.       1.
        9.      10.     324.     295.    1273.    1353.   14843.   13683.       1.
~~~

**`3d_uniaxial_tension.focz`** — applied Z-forces:
~~~bash
      962.   0.33586629E+06
     1003.   0.33586629E+06
     1004.   0.67173257E+06
     1005.   0.67173257E+06
     1006.   0.67173257E+06
~~~

**`3d_uniaxial_tension.disn_1`** — output: 3D nodal displacements (`node_id Ux Uy Uz`, space-separated):
~~~bash
         1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
         2 -0.158000824552E-02  0.000000000000E+00  0.000000000000E+00
         3 -0.599530622114E-04 -0.198353082554E-05  0.000000000000E+00
         4 -0.119662125493E-03 -0.473063658899E-05  0.000000000000E+00
         5 -0.178898410344E-03 -0.864250808527E-05  0.000000000000E+00
         6 -0.237351373110E-03 -0.138088347514E-04  0.000000000000E+00
         7 -0.294744765137E-03 -0.201261481292E-04  0.000000000000E+00
         8 -0.350861002901E-03 -0.273768250648E-04  0.000000000000E+00
         9 -0.405543687855E-03 -0.352738546285E-04  0.000000000000E+00
        10 -0.458697601816E-03 -0.434876069018E-04  0.000000000000E+00
~~~

**`3d_uniaxial_tension_00001.vtk`** — output: 3D VTK snapshot:
~~~bash
# vtk DataFile Version 4.0
X:\PhiPsi_Project\Test_Examples\3d_uniaxial_tension\3d_uniaxial_tension  - results from increment 00001
ASCII
DATASET UNSTRUCTURED_GRID

POINTS 40362 double
    0.000000    0.000000    0.000000
   11.000000    0.000000    0.000000
    0.366667    0.000000    0.000000
    0.733333    0.000000    0.000000
    1.100000    0.000000    0.000000
    ...
CELLS 35152 281216
    8     0     2   119   118  1071  1072  6722  4721
    ...
CELL_TYPES 35152
   12
   12
   ...
POINT_DATA 40362
VECTORS Displacement double
    0.000000E+00 0.000000E+00 0.000000E+00
   -1.580008E-03 0.000000E+00 0.000000E+00
   -5.995306E-05 -1.983531E-06 0.000000E+00
   ...
~~~

<a id="section-Example-D-3d-hf-toughness-dominated"></a>

## Example D — 3D hydraulic fracturing (Test_Examples/3d_hf_toughness_dominated)

A toughness-dominated 3D HF simulation; this example exercises the entire HF output stack.

**`3d_hf_toughness_dominated.kpp`** (key lines):
~~~bash
*Work_Directory
X:\PhiPsi_Project\Test_Examples\3d_hf_toughness_dominated

*Filename
3d_hf_toughness_dominated

*Key_Dimension
3

*Key_Analysis_Type
3                % Hydraulic fracturing.

*Num_Crack
1

*CRACK3D_CIR_COOR_1
...

*INJECTION_T_STAGES_WELLBORES_1_1
*INJECTION_P_STAGES_WELLBORES_1_1
...
~~~

**`3d_hf_toughness_dominated.hftm`** — HF macro/fracture time history (append mode):
~~~bash
    imf   |   ifra   | total_ter|   time   
         1         1         1   0.00000
         1         1         2   0.05000
         1         1         3   0.10000
         1         1         4   0.15000
         1         2         5   0.20000
         1         2         6   0.25000
         2         3         7   0.30000
~~~

**`3d_hf_toughness_dominated.injp`** — injection-pressure history:
~~~bash
  0.000000000000E+00  0.500000000000E+07
  0.500000000000E-01  0.510000000000E+07
  0.100000000000E+00  0.520000000000E+07
  0.150000000000E+00  0.530000000000E+07
  0.200000000000E+00  0.540000000000E+07
~~~

**`3d_hf_toughness_dominated.wbpt`** — wellbore / stage / proppant pressure vs time:
~~~bash
    i_WB   |   i_Stage   | i_Prop   |    Time (s)   |     Pressure (Pa)
         1             1         1       0.00000     0.5000000E+07
         1             1         1       0.05000     0.5100000E+07
         1             1         1       0.10000     0.5200000E+07
         1             1         1       0.15000     0.5300000E+07
         1             1         1       0.20000     0.5400000E+07
~~~

**`3d_hf_toughness_dominated.cms1_1`, `cms2_1`, `cms3_1`** — crack-mesh triangle connectivity (first 5 triangles):
~~~bash
       1
       2
       3
       4
       5
       ...
~~~

(Each value is a 1-based index into the corresponding `.cnox_1`, `.cnoy_1`, `.cnoz_1` arrays.)

**`3d_hf_toughness_dominated.cmap_1`** — crack-mesh aperture (per-crack row):
~~~bash
  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00 ... (50000 values per line, one row per crack)
~~~

**`3d_hf_toughness_dominated.cape_1`** — crack calculation-point aperture:
~~~bash
  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00 ... (2000 values per line, one row per crack)
~~~

**`3d_hf_toughness_dominated.fenx_1`, `feny_1`, `fenz_1`** — fluid-element normal vectors:
~~~bash
  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00  1.000000000000E+00
  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00  1.000000000000E+00
  ...
~~~

---

# Quick-reference summary tables

<a id="section-All-input-files-summary"></a>

## All input-file summary table

| Extension | 2D | 3D | Columns | Required | Reader | Purpose |
|---|---|---|---|---|---|---|
| `.kpp` | Y | Y | text | Yes | `PhiPsi_Read_Input.F90` | Keyword / control input |
| `.node` | Y | Y | 2 / 3 | Yes | `Read_Geo*.F90` | Node coordinates |
| `.elem` | Y | Y | 5 / 9 | Yes | `Read_Geo*.F90` | Element connectivity + material |
| `.boux` | Y | Y | 1 | No | `Read_Geo*.F90` | X-DOF zero constraint |
| `.bouy` | Y | Y | 1 | No | `Read_Geo*.F90` | Y-DOF zero constraint |
| `.bouz` | – | Y | 1 | No | `Read_Geo_3D.F90` | Z-DOF zero constraint |
| `.buxn` | Y | Y | 2 | No | `Read_Geo*.F90` | Non-zero prescribed X disp. |
| `.buyn` | Y | Y | 2 | No | `Read_Geo*.F90` | Non-zero prescribed Y disp. |
| `.buzn` | – | Y | 2 | No | `Read_Geo_3D.F90` | Non-zero prescribed Z disp. |
| `.focx` | Y | Y | 2 | No | `Read_Geo*.F90` | X-direction point load |
| `.focy` | Y | Y | 2 | No | `Read_Geo*.F90` | Y-direction point load |
| `.focz` | – | Y | 2 | No | `Read_Geo_3D.F90` | Z-direction point load |
| `.ivex` / `.ivey` | Y | Y | 2 | No | `Read_Geo*.F90` | Initial velocity X / Y |
| `.ivez` | – | Y | 2 | No | `Read_Geo_3D.F90` | Initial velocity Z |
| `.iacx` / `.iacy` | Y | Y | 2 | No | `Read_Geo*.F90` | Initial acceleration X / Y |
| `.iacz` | – | Y | 2 | No | `Read_Geo_3D.F90` | Initial acceleration Z |
| `.idpx` / `.idpy` | – | Y | 2 | No | `Read_Geo_3D.F90` | Initial displacement X / Y (3D only) |
| `.idpz` | – | Y | 2 | No | `Read_Geo_3D.F90` | Initial displacement Z |
| `.dofx` / `.dofy` | Y | Y | 2 | No | `Read_Geo*.F90` | Multi-point-constraint coupling sets |
| `.fbvl` / `.fbqn` / `.fbiv` | Y | – | 2 | No | `Src_Field_Prob/` | Field-problem boundary value / flux / initial |
| `.bhpc` | Y | – | 2 | Conditional | `Read_Geo*.F90` | Wellbore pressure curve (Key_Analysis_Type 16 / 17) |
| `.eqnl` | – | Y | 1 | Conditional | `Read_Geo_3D.F90` | Earthquake acceleration nodes (EQ_Ac_nodes_list_method = 2) |
| `.boud` / `.focd` | Y | Y | text | No (PPView metadata) | `PPView/Open_Read_Update_kppFile.py` | PPView geometry metadata |

<a id="section-All-output-files-summary"></a>

## All output-file summary table

| Extension | Step suffix | Format | Purpose | Enabling keyword |
|---|---|---|---|---|
| `_<isub5>.vtk` | Yes (zero-padded 5) | ASCII | Main solution snapshot (ParaView) | `*Key_Save_vtk = 1` |
| `_CRACK_<isub5>.vtk` | Yes | ASCII | Crack geometry snapshot | `*Key_Save_vtk = 1` |
| `.disp_<i>` | Yes | ASCII / Binary | Global DOF displacement vector | always (unless `*Key_Save_Nothing = 1`) |
| `.disn_<i>` | Yes | ASCII / Binary | Nodal displacement (`Ux, Uy[, Uz]`) | always |
| `.dien_<i>` | Yes | ASCII / Binary | Enriched-DOF displacement | always |
| `.edei_<i>` | Yes | ASCII | Enriched-DOF index map | always |
| `.disg_<i>` | Yes | ASCII | Gauss-point displacement | always |
| `.dipc_<i>` / `.dinc_<i>` | Yes | ASCII / Binary | Cylindrical-coord displacement | `*Key_CoorSys = 2` |
| `.veln_<i>` / `.acln_<i>` | Yes | ASCII | Per-node velocity / acceleration | dynamic analyses only |
| `.strn_<i>` | Yes | ASCII / Binary | Nodal stress tensor | always |
| `.strg_<i>` | Yes | ASCII / Binary | Gauss-point stress tensor | always |
| `.stnc_<i>` | Yes | ASCII | Nodal stress in cylindrical | `*Key_CoorSys = 2`, 3D |
| `.sttn_<i>` | Yes | ASCII / Binary | Nodal thermal stress | `*Key_Thermal_Stress = 1` |
| `.elss_<i>` | Yes | ASCII | Element 1-3 stress flag | always |
| `.sran_<i>` / `.srac_<i>` | Yes | ASCII / Binary | Nodal strain tensor | always |
| `.gcor_<i>` | Yes | ASCII | Gauss-point coordinates | always |
| `.damg_<i>` | Yes | ASCII | Gauss-point damage factor | cohesive analyses |
| `.elgn_<i>` | Yes | ASCII | Number of Gauss points per element | always |
| `.sifs_<i>` | Yes | ASCII | KI, KII for both tips of every crack | always |
| `.sift_<i>` | Yes | ASCII | SIF time history stack | always |
| `.prst_<i>` | Yes | ASCII | Propagation speed + SIFs at specified tip | `*DY_Save_Propagation_Speed_at_Specified_Tip = 1` |
| `.enee_<i>` | Yes | ASCII | Fracture-energy increment | cohesive analyses |
| `.crax_<i>` / `.cray_<i>` | Yes | ASCII | 2D/3D crack polyline X / Y | always |
| `.craz_<i>` | Yes | ASCII | 3D crack polyline Z | 3D |
| `.crxo_<i>` / `.cryo_<i>` | Yes | ASCII | Original (un-disposed) crack coords | restart runs |
| `.cnox_<i>` / `.cnoy_<i>` / `.cnoz_<i>` | Yes | ASCII | 3D crack-meshed-node coords | 3D |
| `.cms1_<i>` / `.cms2_<i>` / `.cms3_<i>` | Yes | ASCII | Crack-surface triangle connectivity | 3D |
| `.cmso_<i>` | Yes | ASCII | Crack-mesh outline | 3D |
| `.cmap_<i>` | Yes | ASCII | Crack-mesh aperture (3D) | 3D |
| `.cape_<i>` | Yes | ASCII | Crack calc-point aperture | HF / cohesive |
| `.capf_<i>` | Yes | ASCII | Aperture at fluid elements (3D HF) | 3D HF |
| `.ctap_<i>` | Yes | ASCII | Crack tangential opening | cohesive / contact |
| `.cvpx_<i>` … `.cvzz_<i>` | Yes | ASCII | Crack vertex coords + 9-axis tangents | 3D |
| `.apex_<i>` / `.apey_<i>` | Yes | ASCII | Calc-point coordinates | 2D HF |
| `.cori_<i>` | Yes | ASCII | Calc-point orientation angles | 2D HF |
| `.cpre_<i>` / `.cvel_<i>` / `.cqua_<i>` | Yes | ASCII | Calc-point pressure / flow velocity / quantity | HF |
| `.cohx_<i>` / `.cohy_<i>` | Yes | ASCII | Cohesive tractions X / Y | cohesive |
| `.crrd_<i>` | Yes | ASCII | Crack radius (radial cracks) | `*Key_Save_Crack_Radius = 1` |
| `.ennd_<i>` / `.enns_<i>` / `.ennh_<i>` / `.ennj_<i>` | Yes | ASCII | Enriched-node type per crack / cross / hole / inclusion | always |
| `.elty_<i>` / `.elts_<i>` / `.elth_<i>` / `.eltj_<i>` | Yes | ASCII | Element type per crack / cross / hole / inclusion | always |
| `.posi_<i>` / `.poss_<i>` / `.posh_<i>` / `.posj_<i>` | Yes | ASCII | Signed-distance sign per node per feature | always |
| `.njel_<i>` / `.nods_<i>` | Yes | ASCII | Crack-node / host-element associations | always |
| `.celt_<i>` / `.celv_<i>` / `.celj_<i>` / `.celc_<i>` | Yes | ASCII | Tip / vertex / junction / intersection coords | always |
| `.ctty_<i>` | Yes | ASCII | Crack-tip type | always |
| `.blab_<i>` | Yes | ASCII | Crack baseline vectors | always |
| `.blvx_<i>` / `.blvy_<i>` / `.blvz_<i>` | Yes | ASCII | Baseline direction-vector components | always |
| `.tere_<i>` | Yes | ASCII | Tip-enriched node element number | always |
| `.cndx_<i>` / `.cndy_<i>` / `.cndz_<i>` | Yes | ASCII | Crack-node S1 vector (CFCP-2) | `*CFCP = 2` |
| `.hlcr` | No | ASCII | Circular hole coordinates | holes enabled |
| `.ehcr` | No | ASCII | Elliptical hole parameters | holes enabled |
| `.jzcr` | No | ASCII | Circular inclusion coordinates | inclusions enabled |
| `.jzpx` / `.jzpy` | No | ASCII | Polygonal inclusion vertices | inclusions enabled |
| `.cscr` | No | ASCII | Cross-interface coordinates | cross enabled |
| `.nfcx` / `.nfcy` / `.nfcz` | No | ASCII | Natural-fracture vertex coordinates (3D) | natural fractures |
| `.ncrx` / `.ncry` | No | ASCII | Natural-crack endpoint coords (2D) | natural fractures |
| `.ener` | Append | ASCII | Energy balance per step | always |
| `.edye` | Append | ASCII | Explicit-dynamics energy history | explicit dynamics |
| `.edtm` | Append | ASCII | Explicit-dynamics time history | explicit dynamics |
| `.edap` | Append | ASCII | Explicit-dynamics aperture history (3D) | explicit dynamics |
| `.dcrl` / `.dprl` | Append | ASCII | Crack total / propagation length history | dynamic 2D |
| `.idtm` | Append | ASCII | Implicit-dynamics time history | implicit dynamics |
| `.iite` / `.ilth` / `.ipre` | Append | ASCII | iFrac iterations / length / pressure | HF |
| `.lpmf_<i>` / `.lpmx_<i>` / `.lpma_<i>` | Yes | ASCII | Lumped mass per node (FEM / enriched / both) | dynamic |
| `.hftm` | Append | ASCII | HF macro/fracture time log | HF |
| `.injp` | Append | ASCII | Injection-pressure history | HF |
| `.wbfp` | No | ASCII | Wellbore fracturing path | HF |
| `.wbpt` | Append | ASCII | Wellbore / stage / proppant pressure vs time | HF |
| `.ihft` | No | ASCII | HF time per fracture step | HF |
| `.icpt` / `.icpt_LS` / `.icpt_CT` | No | ASCII | CPU-time per step (total / linear-solver / contact) | HF |
| `.ccon_<i>` / `.pokf_<i>` / `.cond_<i>` / `.wpnp_<i>` / `.epcr_<i>` | Yes | ASCII | Proppant transport outputs | `*Key_Propp_Trans = 1` |
| `.Saved_Filename` | Append | ASCII | Flow rate at specific point | HF |
| `.elcs_<i>` / `.elco_<i>` | Yes | ASCII | Element contact / cohesive state | contact / cohesive |
| `.kiel_<i>` | Yes | ASCII | Killed (broken) element IDs | element deactivation |
| `.fxdf_<i>` / `.fydf_<i>` / `.fzdf_<i>` | Yes | ASCII | DOF internal force vectors | always |
| `.fxsl_<iSL>_<i>` / `.fysl_<iSL>_<i>` / `.fzsl_<iSL>_<i>` | Yes | ASCII | Surface-load internal force vectors | 3D surface loads |
| `.skxf` | No | Binary | XFEM element stiffness matrices | always (debug) |
| `.csrn_<i>` / `.csra_<i>` / `.csrj_<i>` / `.csri_<i>` | Yes | ASCII | CSR stiffness matrix | always (debug) |
| `.fenx_<i>` … `.fezz_<i>` | Yes | ASCII | Fluid-element normals / up-low / local axes (3D HF) | 3D HF |
| `.fnnx_<i>` / `.fnny_<i>` / `.fnnz_<i>` | Yes | ASCII | Fluid-node normals | 3D HF |
| `.fnux_<i>` / `.fnuy_<i>` / `.fnuz_<i>` / `.fnlx_<i>` / `.fnly_<i>` / `.fnlz_<i>` | Yes | ASCII | Up-face / low-face crack-node displacements | 3D HF |
| `.cpfn_<i>` / `.cpno_<i>` | Yes | ASCII | Fluid-node count / numbers per crack | 3D HF |
| `.ccpx_<i>` / `.ccpy_<i>` / `.ccpz_<i>` | Yes | ASCII | Fluid-node coordinates | 3D HF |
| `.cnlx_<i>` / `.cnly_<i>` / `.cnlz_<i>` | Yes | ASCII | Crack-node local coords | 3D HF |
| `.cmse_<i>` | Yes | ASCII | Element containing each crack node | 3D HF |
| `.fraz` | No | ASCII | Fracture-zone bounding box | always |
| `.seed` | No | ASCII | Random-number seed | stochastic analyses |
| `.post` | No | ASCII | MATLAB post-processor header | always |
| `.fdcu` / `.fccu` | No | ASCII | Force-displacement / Force-COD curves | displacement control |
| `.rbco_<i>` | Yes | ASCII | Rigid-ball (circle) coordinates | rigid-ball contact |
| `.fdvl_<i>` | Yes | ASCII | Per-node fluid velocity | HF |
| `.ecfv_<i>` | Yes | ASCII | Element-centroid field value | HF |
| `.sccx_<i>` / `.sccy_<i>` / `.scdx_<i>` / `.scdy_<i>` | Yes | ASCII | Stress-corrosion-cracking coords / displacements | SCC analyses |
| `PhiPsi_Console_Window.log` | Append | text | Console-window capture | always |
| `current_folder.dat` | No | text | Last-used `Work_Directory` | always |

---

# End of manual

For keyword reference, see [PhiPsi Keywords Manual](http://www.phipsi.top/phipsi_keywords_manual.html). For the project homepage, see <http://phipsi.top>.

If you find a missing file extension or a format that disagrees with the PhiPsi data, please report it to the author (see email at top).