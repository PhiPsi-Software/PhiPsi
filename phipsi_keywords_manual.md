# PhiPsi Keywords Manual
+ This manual simply describes the usage of each keyword available in PhiPsi.
+ Author: Fang Shi
+ Website: http://phipsi.top
+ Email: shifang@hyit.edu.cn / shifang@ustc.edu.cn
+ Updated on 2026-02-28.

---

Click [here](http://phipsi.top/source/PhiPsi_Win64_Latest.zip) to download the latest version of PhiPsi.

Click [here](https://sourceforge.net/projects/phipsi/files/) to download the latest version of PPView.

---

# Contents of this keywords manual

[Keywords rules](#section-Keywords-rules) 

[Basic control keywords](#section-Basic-control-keywords)

[Linear solvers](#section-Linear-solvers)

[Define initial cracks, voids (holes) and inclusions](#section-Define-initial-cracks-voids-holes-and-inclusions)

[Definition of material parameters](#section-Definition-of-material-parameters)

[Hydraulic fracturing simulation](#section-Hydraulic-fracturing-simulation)

[Nonlinear analysis](#section-Nonlinear-analysis)

[Cohesive crack](#section-Cohesive-crack) 

[Dynamic analysis](#section-Dynamic-analysis)

[Coupling of degrees of freedom](#section-Coupling-of-degrees-of-freedom)

[Deactivation elements](#section-Deactivation-elements)

[Thermal load and prestress](#section-Thermal-load-prestress)

[Surface load](#section-Surface-load)

[Control of the program](#section-Control-of-the-program)

[Keywords file examples](#section-Keywords-file-examples)

<a id="section-Keywords-rules"></a>

# Keywords rules

+ Keyword lines are in no particular order.
  
+ Lines starting with "%" are comments.

+ The keyword line starts with "*", followed by a corresponding data line (multiple data are separated by commas ","). For example:
  
    ~~~bash
    *num_Crack
    2
    *InSitu_S1_nv_3D
    1.0,0.0,0.0
    ~~~

+ Keywords are not case sensitive.

+ Support parameter definition. For example:
  
    ~~~bash
    num_of_crack = 2
    *num_Crack
    num_of_crack
    ~~~
  
+ Support basic four arithmetic operations ($+$, $-$, $*$, $()$, and $/$). For example:
  
    ~~~bash
    *num_Crack
    1+1

    *INI_CRACK_PRESSURE_1
    (15.0e6-10.0e6)*(3.0-1.0)

    Time_Hour=1.0
    *INJECTION_T_STAGES_WELLBORES_1_1
    Time_hour*60*60
    ~~~

<a id="section-Basic-control-keywords"></a>

# Basic control keywords

+  ***<font color=green>Work_Dirctory</font>** - File location of the input files of PhiPsi, i.e., the work directory. 
    
    Keyword example:
    ~~~bash
    *Work_Directory
    X:\PhiPsi_Work\FEM
    ~~~
+ ***<font color=green>Filename</font>** - Filename of the input files of PhiPsi (all input files have the same filename).

    Keyword example:    
    ~~~bash
    *Filename
    FEM
    ~~~
+ ***<font color=green>Key_Unit_System</font>** - Unit system. 
    
    = 1, international system of units (default); 

    = 2, mm-ton-s.
    
    Keyword example:
    ~~~bash
    *Key_Unit_System
    1     % SI.
    ~~~           
+ ***<font color=green>Key_Dimension</font>** - The dimension of the problem. 
    
    = 2, two-dimension (Plane stress or plane strain)
    
    = 3, three-dimension

    Keyword example:
    ~~~bash
    *Key_Dimension
    2     % 2D problem.
    ~~~ 
+ ***<font color=green>Key_Type_2D</font>** - Define plane stress or plane strain for 2D problems.
  
    = 1, plane stress;

    = 2, plane strain (default). 
    
    Keyword example:
    ~~~bash
    *Key_Type_2D
    2     % plane strain.
    ~~~ 
+ ***<font color=green>Key_Analysis_Type</font>** - Analysis type.

    = 1, quasi-static (default); 

    = 2, implicit dynamic;

    = 3, hydraulic fracturing simulation;

    = 4, hydraulic fracturing simulation with slick water;

    = 7, nonlinear problem (such as plastic deformation problem); 

    = 15, field problem (such as heat transfer problem). 

    Keyword example:
    ~~~bash
    *Key_Analysis_Type
    1     % quasi-static.
    ~~~ 
+ ***<font color=green>Key_SIFs_Method</font>** - Method to calculate the stress intensity factors (SIFs).

    = 1, displacement interpolation method (DIM; default); 

    = 2, the interaction integral method.

    Keyword example:
    ~~~bash
    *Key_SIFs_Method
    2     % the interaction integral method.
    ~~~ 
+ ***<font color=green>Key_SIFs_DIM_Points</font>** - Number of points to calculate SIFs using DIM method (default to 2).

    Keyword example:
    ~~~bash
    *Key_SIFs_DIM_Points
    2  
    ~~~     
+ ***<font color=green>Key_SIFs_DIM_Method</font>** - Model assumption of the DIM method.

    = 1, plane stress (default); 

    = 2, plane strain.

    Keyword example:
    ~~~bash
    *Key_SIFs_DIM_Method
    2     
    ~~~     
+ ***<font color=green>Key_Print_SIFs_to_Screen</font>** - Print SIFs to screen / terminal.

    = 0, do not print SIFs to screen / terminal (default); 

    = 1, print SIFs to screen / terminal.

    Keyword example:
    ~~~bash
    *Key_Print_SIFs_to_Screen
    1     
    ~~~         
+ ***<font color=green>Key_Contact</font>** - Consider contact between crack surfaces or not.
  
    = 0, no (default); 
    
    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Contact
    1     % Active contact between crack surfaces.
    ~~~ 
+ ***<font color=green>Max_Contact_Iter</font>** - The maximum number of contact iterations (default to 10).

    Keyword example:
    ~~~bash
    *Max_Contact_Iter
    10    
    ~~~     
+ ***<font color=green>fric_mu_Cont</font>** - Friction coefficient of crack surfaces (default to 0.3). 

    Keyword example:
    ~~~bash
    *fric_mu_Cont
    0.4
    ~~~ 
+ ***<font color=green>kn_Cont_Penalty</font>** - Normal penalty stiffness of crack surfaces (default to 1.0e13). 

    Keyword example:
    ~~~bash
    *kn_Cont_Penalty
    1.0e12
    ~~~ 
+ ***<font color=green>kt_Cont_Penalty</font>** - Tangential penalty stiffness of crack surfaces (default to 1.0e13). 

    Keyword example:
    ~~~bash
    *kt_Cont_Penalty
    1.0e12
    ~~~ 
+ ***<font color=green>Conve_Tol_Penalty</font>** - Convergence tolerance for the contact iteration (default to 1.0e-5). 

    Keyword example:
    ~~~bash
    *Conve_Tol_Penalty
    1.0e-4    
    ~~~ 
+ ***<font color=green>Num_Substeps</font>** - Number of steps need to be performed.

    Keyword example:
    ~~~bash
    *Num_Substeps
    10     
    ~~~ 
+ ***<font color=green>CFCP</font>** - Criterion of crack propagation.

    = 1, maximum tensile circumferential stress criterion (default); 

    = 2, maximum principal tensile stress criterion;
    
    = 3, Schollmann’s criterion;

    = 4, maximum shear stress criterion.

    Keyword example:
    ~~~bash
    *CFCP
    2      % Maximum principal tensile stress criterion.
    ~~~ 
+ ***<font color=green>Boundary_Cracks</font>** - Mark 3D cracks as edge (boundary) cracks.
  
    = 0, internal crack; 

    = 1, edge (boundary) cracks.

    Keyword example:
    ~~~bash
    *Boundary_Cracks
    1,1,0    % Mark cracks 1 and 2 as edge (boundary) cracks.  
    ~~~     
+ ***<font color=green>Key_CFCP_3_Type</font>** - Crack growth behavior for the Schollmann’s criterion.
  
    = 1, crack grows only when the equivalent stress intensity factor is greater than the fracture toughness $K_{Ic}$ (default); 

    = 2, crack grows when the equivalent stress intensity factor is greater than zero.

    Keyword example:
    ~~~bash
    *Key_CFCP_3_Type
    2     
    ~~~ 
+ ***<font color=green>Key_Ave_Stress</font>** - Calculation method of weighted average stress.

    = 1, Shi's Formation, See p25 of Study on the Cracking Process of Rock Using the Extended Finite Element Method; 

    = 2, Ordinary Formation, See p498 of Extended Finite Element Method: Theory and Applications.

    Keyword example:
    ~~~bash
    *Key_Ave_Stress
    1     
    ~~~     
+ ***<font color=green>FACTOR_AVE_R</font>** - Radius factor for calculating weighted average stress (default to 0.7).
  
    Keyword example:
    ~~~bash
    *FACTOR_AVE_R
    0.7     
    ~~~    
+ ***<font color=green>A_AVE_SHI</font>** - Parameter a for calculating weighted average stress when Key_Ave_Stress = 1 (1-100, the larger it is, the closer it gets to the average; default to 5).
  
    Keyword example:
    ~~~bash
    *A_AVE_SHI
    5  
    ~~~          
+ ***<font color=green>Key_TipEnrich</font>** - Crack tip enrichment type.

    = 0, no crack tip enrichment, the crack tip is automatically adjusted to the element boundary; 

    = 1, standard crack tip enrichment (4 items, $F_1$ to $F_4$; default).

    = 2, keep only the first item ($F_1$, this is recommended for dynamic analysis); 

    = 3, not available yet;

    = 4, Cohesive crack tip_cohesive crack (only 1 enrichment function item).

    Keyword example:
    ~~~bash
    *Key_TipEnrich
    1     
    ~~~     
+ ***<font color=green>Key_Fd_TipEnrich</font>** - Crack tip enrichment type for field problems.

    = 0, no crack tip enrichment, the crack tip is automatically adjusted to the element boundary; 

    = 1, strong discontinuity crack tip enrichment (enrichment function is $sqrt(r)*sin(θ/2)$; adiabatic crack for heat conduction problem; default);

    = 2, weak discontinuity crack tip enrichment (enrichment function is $sqrt(r)*cos(θ/2)$; isothermal crack for heat conduction problem).

    Keyword example:
    ~~~bash
    *Key_Fd_TipEnrich
    1
    ~~~      
+ ***<font color=green>Key_InPlane_Growth</font>** - 3D Plannar or nonplanar crack propagation.

    = 0, nonplanar crack propagation (default); 

    = 1, planar crack propagation.

    Keyword example:
    ~~~bash
    *Key_InPlane_Growth
    1 
    ~~~    
+ ***<font color=green>Key_Scheme_Signed_Dis_InPlane</font>** - Scheme to calculate the signed distance for 3D plannar crack.

    = 1, slow and robuts (default); 

    = 2, fast but not robuts.
    
    = 3, more fast but not robuts.
    Keyword example:
    ~~~bash
    % Fast scheme to calculte the signed distance to plannar crack.
    *Key_Scheme_Signed_Dis_InPlane
    2 
    ~~~     
+ ***<font color=green>Key_Allow_3D_Outside_Crack</font>** - Allow the 3D cracks to be outside the model. If Boundary_Cracks is set to 1 for any 3D crack, the program automatically sets Key_Allow_3D_Outside_Crack to 1.

    = 0, does not allow the cracks to be outside the model (default); 

    = 1, allow the cracks to be outside the model.

    Keyword example:
    ~~~bash
    *Key_Allow_3D_Outside_Crack
    1 
    ~~~      
 + ***<font color=green>Key_Stop_Outside_Crack</font>** - Do not allow cracks to grow outside the model.

    = 0, allow cracks to grow outside the model (default); 

    = 1, do not allow cracks to grow outside the model.

    Keyword example:
    ~~~bash
    *Key_Stop_Outside_Crack
    1 
    ~~~     
+ ***<font color=green>Key_Denoise_Vertex_Value</font>** - Denoise vertex values for 3D cracks.

    = 0, no (default); 

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Denoise_Vertex_Value
    1 
    ~~~      
+ ***<font color=green>Key_Smooth_Vertex_Value</font>** - Smooth vertex values for 3D cracks.

    = 0, no (default); 

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Smooth_Vertex_Value
    1 
    ~~~     
+ ***<font color=green>Key_Smooth_GF_Value</font>** - Smooth the propagation for 3D cracks.

    = 0, no (default); 

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Smooth_GF_Value
    1 
    ~~~  
+ ***<font color=green>Key_Smooth_Front</font>** - 3D crack propagation front vertex smoothing scheme.

    = 0, do not smooth crack front (default);

    = 1, Smooth the crack front edge according to the straight line equation;

    = 2, According to the equation of the circle;

    = 3, CSAPS cubic spline fitting;

    = 4, B-spline curve fitting;

    = 5, According to the equation of ellipse;

    = 6, Taubin smoothing algorithm.

    Keyword example:
    ~~~bash
    *Key_Smooth_Front
    6 
    ~~~  
+ ***<font color=green>Key_Force_Control</font>** - The method to control the value of the applied force for each step.

    = 1, the force is applied all at once (default); 

    = 2, the applied force increases linearly in each step;

    = 3, not available yet;

    = 4, special scheme for cohesive crack;

    = 5, make sure that only one crack grows in each step. 

    Keyword example:
    ~~~bash
    *Key_Force_Control
    1      
    ~~~       
+ ***<font color=green>Key_Initiation</font>** - Allow the formation of new cracks.

    = 0, no (default);

    = 1, yes.

    Keyword example:
    ~~~bash
    % New cracks are allowed to be generated.
    *Key_Initiation
    1      
    ~~~ 
+ ***<font color=green>Key_Ini_Rule</font>** - Ruls for the formation of new cracks.

    = 1, maximum tensile stress criterion (default);

    = 2, maximum shear stress criterion;

    = 3, concrete.

    Keyword example:
    ~~~bash
    % New cracks are allowed to be generated.
    *Key_Ini_Rule
    1      
    ~~~     
+ ***<font color=green>Key_Ini_Cr_3D_Type</font>** - Shape of newly formed cracks.

    = 1, circular (default);

    = 2, rectangle.

    Keyword example:
    ~~~bash
    *Key_Ini_Cr_3D_Type
    1      
    ~~~   
+ ***<font color=green>Size_Ini_Crack</font>** - Size of newly formed cracks. For a 3D circular crack, Size_Ini_Crack represnts the diameter.

    Keyword example:
    ~~~bash
    *Size_Ini_Crack
    0.5     
    ~~~   
+ ***<font color=green>Key_Max_Num_Initiation_Cracks</font>** - Maximum number of cracks allowed to be newly formed (default to 1).

    Keyword example:
    ~~~bash
    *Key_Max_Num_Initiation_Cracks
    10     
    ~~~   
+ ***<font color=green>Key_Propagation</font>** - Propagation of cracks.

    = 0, cracks are not allowed to propagate;

    = 1, cracks are allowed to propagate (default).

    Keyword example:
    ~~~bash
    % Cracks are allowed to propagate.
    *Key_Propagation
    1      
    ~~~     
+ ***<font color=green>MAT_ALLOW_CRACK_INITIATION_1 - MAT_ALLOW_CRACK_INITIATION_20</font>** - Allow the formation of new cracks in material type $n$.

    = 0, do not allow;

    = 1, allow (default).

    Keyword example:
    ~~~bash
    % Do not allow the formation of new cracks in material type 2.
    *MAT_ALLOW_CRACK_INITIATION_2
    0     
    ~~~   
+ ***<font color=green>Cracks_Allow_Propa</font>** - Allow cracks 1 to $n$ to propagate.

    = 0, does not allow cracks to propagate; 

    = 1, allow cracks to propagate (default).

    Keyword example:
    ~~~bash
    *Cracks_Allow_Propa
    1,1,1,1,1,0,0,0,0,0     % Allow cracks 1 to 5 to propagate, and does not allow the propagation of cracks 6 to 10.
    ~~~      
+ ***<font color=green>Key_CS_Natural_Crack</font>** - Using penalty function method to simulate contact of natural cracks.

    = 0, no (default); 

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_CS_Natural_Crack
    1
    ~~~     
+ ***<font color=green>Key_Check_and_Adjust_Cracks_3D</font>** - Check and adjust initial 3D cracks.

    = 0, no (default); 

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Check_and_Adjust_Cracks_3D
    1
    ~~~    
+ ***<font color=green>Factor_Propagation</font>** - Factor of propagation length of cracks, i.e., propagation length $\Delta l=*Factor_Propagation \times l_c$, where $l_c$ represents the average size of enriched elements. 

    Keyword example:
    ~~~bash
    *Factor_Propagation
    1.5      
    ~~~     
+ ***<font color=green>Propagation_Length</font>** - Propagation length of cracks during a propagation step (Unit: $m$). 

    Keyword example:
    ~~~bash
    *Propagation_Length
    0.010   
    ~~~        
+ ***<font color=green>Key_Fracture_Zone</font>** - Define propagation zone, i.e., cracks are allowed to propagate just inside the propagation zone.

    = 0, no propagation zone (default);

    = 1, define a propagation zone.

    Keyword example:
    ~~~bash
    *Key_Fracture_Zone
    1      
    ~~~     
+ ***<font color=green>Frac_Zone_MinX</font>** - Minimum value of range $x$ of the fracture zone.

    Keyword example:
    ~~~bash
    *Frac_Zone_MinX
    4.0e-3    
    ~~~       
+ ***<font color=green>Frac_Zone_MinY</font>** - Minimum value of range $y$ of the fracture zone.

    Keyword example:
    ~~~bash
    *Frac_Zone_MinY
    4.0e-3    
    ~~~    
+ ***<font color=green>Frac_Zone_MinZ</font>** - Minimum value of range $z$ of the fracture zone.

    Keyword example:
    ~~~bash
    *Frac_Zone_MinZ
    4.0e-3    
    ~~~        
+ ***<font color=green>Frac_Zone_MaxX</font>** - Maximum value of range $x$ of the fracture zone.

    Keyword example:
    ~~~bash
    *Frac_Zone_MaxX
    50.0e-3   
    ~~~       
+ ***<font color=green>Frac_Zone_MaxY</font>** - Maximum value of range $y$ of the fracture zone.

    Keyword example:
    ~~~bash
    *Frac_Zone_MaxY
    50.0e-3   
    ~~~    
+ ***<font color=green>Frac_Zone_MaxZ</font>** - Maximum value of range $z$ of the fracture zone.

    Keyword example:
    ~~~bash
    *Frac_Zone_MaxZ
    50.0e-3   
    ~~~  
+ ***<font color=green>Key_Ini_Crack_Zone</font>** - Define a zone for crack initiation, i.e., initial cracks are allowed to be created just inside the given zone.

    = 0, no (default);

    = 1, define a crack initiation zone.

    Keyword example:
    ~~~bash
    *Key_Ini_Crack_Zone
    1      
    ~~~     
+ ***<font color=green>Ini_Crack_Zone_MinX</font>** - Minimum value of range $x$ of the crack initiation zone.

    Keyword example:
    ~~~bash
    *Ini_Crack_Zone_MinX
    4.0e-3    
    ~~~       
+ ***<font color=green>Ini_Crack_Zone_MinY</font>** - Minimum value of range $y$ of the crack initiation zone.

    Keyword example:
    ~~~bash
    *Ini_Crack_Zone_MinY
    4.0e-3    
    ~~~    
+ ***<font color=green>Ini_Crack_Zone_MaxX</font>** - Maximum value of range $x$ of the crack initiation zone.

    Keyword example:
    ~~~bash
    *Ini_Crack_Zone_MaxX
    4.0e-3    
    ~~~        
+ ***<font color=green>Ini_Crack_Zone_MaxY</font>** - Maximum value of range $y$ of the crack initiation zone.

    Keyword example:
    ~~~bash
    *Ini_Crack_Zone_MaxY
    50.0e-3   
    ~~~       
+ ***<font color=green>Key_InSitu_Strategy</font>** - In-situ stress.

    = 0, does not consider initial stress issues;

    = 1, the displacement field generated by in-situ stress is processed through linear superposition;

    = 2, use the method proposed in the book 'The Finite Element Method: its Basis and Fundamentals';

    = 3, shield the applied in-situ stress (only for hydraulic fracturing analysis);

    = 4, specify the initial stress or read from the initial stress file.

    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4   
    ~~~  
+ ***<font color=green>Key_Read_Initial_Node_Stress_File</font>** - Read In-situ stress from the initial stress file (.istn file) when Key_InSitu_Strategy = 4. *.istn file is the initial stress file of all nodes. The number of rows of this file is the same as the total number of nodes. There are 6 data per row: $S_{xx}$, $S_{yy}$, $S_{zz}$, $S_{xy}$, $S_{yz}$, $S_{xz}$ and the value of stress is positive in tension.

    = 0, define In-situ stress manualy (default);

    = 1, read In-situ stress from the initial stress file (*.istn file).

    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4   
    *Key_Read_Initial_Node_Stress_File
    1
    ~~~   
+ ***<font color=green>InSitu_S1_3D</font>** - In-situ stress $S_1$ when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *InSitu_S1_3D
    10.0e6    
    ~~~   
+ ***<font color=green>InSitu_S1_nv_3D</font>** - Normal vector of In-situ stress $S_1$ when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *InSitu_S1_3D
    10.0e6    
    *InSitu_S1_nv_3D
    1.0,0.0,0.0
    ~~~            
+ ***<font color=green>InSitu_S2_3D</font>** - In-situ stress $S_2$ when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *InSitu_S2_3D
    20.0e6    
    ~~~   
+ ***<font color=green>InSitu_S2_nv_3D</font>** - Normal vector of In-situ stress $S_2$ when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *InSitu_S2_3D
    20.0e6 
    *InSitu_S2_nv_3D
    0.0,1.0,0.0
    ~~~         
+ ***<font color=green>InSitu_S3_3D</font>** - In-situ stress $S_3$ when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *InSitu_S3_3D
    5.0e6    
    ~~~   
+ ***<font color=green>InSitu_S3_nv_3D</font>** - Normal vector of In-situ stress $S_3$ when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *InSitu_S3_3D
    5.0e6  
    *InSitu_S3_nv_3D
    0.0,0.0,1.0
    ~~~         
+ ***<font color=green>Key_Nonuniform_InSitu_X_with_Z</font>** - Apply Nonuniform InSitu stress in X direction (changing with Z) when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_X_with_Z
    1
    ~~~          
+ ***<font color=green>InSitu_Sx_3D_Seg_Strs_X_with_Z</font>** - Sx stresses when Key_Nonuniform_InSitu_X_with_Z=1 and Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    % In this example, three different stress zones (with different stress values of 15, 10, and 5 MPa, respectively) are defined.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_X_with_Z
    1
    *InSitu_Sx_3D_Seg_Strs_X_with_Z
    15.0e6,10.0e6,5.0e6
    ~~~        
+ ***<font color=green>InSitu_Sx_3D_Seg_Loca_X_with_Z</font>** - Location of Sx stresses when Key_Nonuniform_InSitu_X_with_Z=1 and Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    % In this example, along the z-direction, the stress-x in the zone (0.0, 110.0) is 15.0 MPa, the stress-x in the zone (110.0, 140.0) is 10.0 MPa, and the stress-x in the zone (140.0, 250.0) is 5.0 MPa.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_X_with_Z
    1
    *InSitu_Sx_3D_Seg_Strs_X_with_Z
    15.0e6,10.0e6,5.0e6
    *InSitu_Sx_3D_Seg_Loca_X_with_Z
    0.0,110.0,140.0,250.0   
    ~~~    
+ ***<font color=green>Key_Nonuniform_InSitu_X_with_Y</font>** - Apply nonuniform InSitu stress in X direction (changing with Y) when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_X_with_Z
    1
    ~~~          
+ ***<font color=green>InSitu_Sx_3D_Seg_Strs_X_with_Y</font>** - Sx stresses when Key_Nonuniform_InSitu_X_with_Y=1 and Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    % In this example, three different stress zones (with different stress values of 15, 10, and 15 MPa, respectively) are defined.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_X_with_Z
    1
    *InSitu_Sx_3D_Seg_Strs_X_with_Y
    15.0e6,10.0e6,15.0e6
    ~~~        
+ ***<font color=green>InSitu_Sx_3D_Seg_Loca_X_with_Y</font>** - Locations of Sx stresses when Key_Nonuniform_InSitu_X_with_Y=1 and Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    % In this example, along the y-direction, the stress-x in the zone (0.0, 110.0) is 15.0 MPa, the stress-x in the zone (110.0, 140.0) is 10.0 MPa, and the stress-x in the zone (140.0, 250.0) is 15.0 MPa.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_X_with_Y
    1
    *InSitu_Sx_3D_Seg_Strs_X_with_Y
    15.0e6,10.0e6,15.0e6
    *InSitu_Sx_3D_Seg_Loca_X_with_Y
    0.0,110.0,140.0,250.0   
    ~~~   
+ ***<font color=green>Key_Nonuniform_InSitu_Y_with_Z</font>** - Apply nonuniform InSitu stress in Y direction (changing with Z) when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Y_with_Z
    1
    ~~~  
+ ***<font color=green>InSitu_Sy_3D_Seg_Strs_Y_with_Z</font>** - Sy stresses when Key_Nonuniform_InSitu_Y_with_Z=1 and Key_InSitu_Strategy = 4.   
    Keyword example:
    ~~~bash
    % In this example, three different stress zones (with different stress values of 15, 10, and 15 MPa, respectively) are defined.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Y_with_Z
    1
    *InSitu_Sy_3D_Seg_Strs_Y_with_Z
    15.0e6,10.0e6,15.0e6
    ~~~   
+ ***<font color=green>InSitu_Sy_3D_Seg_Loca_Y_with_Z</font>** - Locations of Sy stresses when Key_Nonuniform_InSitu_Y_with_Z=1 and Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    % In this example, along the z-direction, the stress-y in the zone (0.0, 110.0) is 15.0 MPa, the stress-y in the zone (110.0, 140.0) is 10.0 MPa, and the stress-y in the zone (140.0, 250.0) is 15.0 MPa.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Y_with_Z
    1
    *InSitu_Sy_3D_Seg_Strs_Y_with_Z
    15.0e6,10.0e6,15.0e6
    *InSitu_Sy_3D_Seg_Loca_Y_with_Z
    0.0,110.0,140.0,250.0   
    ~~~   
+ ***<font color=green>Key_Nonuniform_InSitu_Y_with_X</font>** - Apply nonuniform InSitu stress in Y direction (changing with X) when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Y_with_X
    1
    ~~~  
+ ***<font color=green>InSitu_Sy_3D_Seg_Strs_Y_with_X</font>** - Sy stresses when Key_Nonuniform_InSitu_Y_with_X=1 and Key_InSitu_Strategy = 4.     
    Keyword example:
    ~~~bash
    % In this example, three different stress zones (with different stress values of 15, 10, and 15 MPa, respectively) are defined.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Y_with_X
    1
    *InSitu_Sy_3D_Seg_Strs_Y_with_X
    15.0e6,10.0e6,15.0e6
    ~~~   
+ ***<font color=green>InSitu_Sy_3D_Seg_Loca_Y_with_X</font>** - Locations of Sy stresses when Key_Nonuniform_InSitu_Y_with_X=1 and Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    % In this example, along the x-direction, the stress-y in the zone (0.0, 110.0) is 15.0 MPa, the stress-y in the zone (110.0, 140.0) is 10.0 MPa, and the stress-y in the zone (140.0, 250.0) is 15.0 MPa.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Y_with_X
    1
    *InSitu_Sy_3D_Seg_Strs_Y_with_X
    15.0e6,10.0e6,15.0e6
    *InSitu_Sy_3D_Seg_Loca_Y_with_X
    0.0,110.0,140.0,250.0   
    ~~~   
+ ***<font color=green>Key_Nonuniform_InSitu_Z_with_X</font>** - Apply nonuniform InSitu stress in Z direction (changing with X) when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Z_with_X
    1
    ~~~  
+ ***<font color=green>InSitu_Sz_3D_Seg_Strs_Z_with_X</font>** - Sz stresses when Key_Nonuniform_InSitu_Z_with_X=1 and Key_InSitu_Strategy = 4.     
    Keyword example:
    ~~~bash
    % In this example, three different stress zones (with different stress values of 15, 10, and 15 MPa, respectively) are defined.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Z_with_X
    1
    *InSitu_Sz_3D_Seg_Strs_Z_with_X
    15.0e6,10.0e6,15.0e6
    ~~~   
+ ***<font color=green>InSitu_Sz_3D_Seg_Loca_Z_with_X</font>** - Locations of Sz stresses when Key_Nonuniform_InSitu_Z_with_X=1 and Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    % In this example, along the x-direction, the stress-z in the zone (0.0, 110.0) is 15.0 MPa, the stress-y in the zone (110.0, 140.0) is 10.0 MPa, and the stress-y in the zone (140.0, 250.0) is 15.0 MPa.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Z_with_X
    1
    *InSitu_Sz_3D_Seg_Strs_Z_with_X
    15.0e6,10.0e6,15.0e6
    *InSitu_Sz_3D_Seg_Loca_Z_with_X
    0.0,110.0,140.0,250.0   
    ~~~   
+ ***<font color=green>Key_Nonuniform_InSitu_Z_with_Y</font>** - Apply nonuniform InSitu stress in Z direction (changing with Y) when Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Z_with_Y
    1
    ~~~  
+ ***<font color=green>InSitu_Sz_3D_Seg_Strs_Z_with_Y</font>** - Sz stresses when Key_Nonuniform_InSitu_Z_with_Y=1 and Key_InSitu_Strategy = 4.     
    Keyword example:
    ~~~bash
    % In this example, three different stress zones (with different stress values of 15, 10, and 15 MPa, respectively) are defined.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Z_with_Y
    1
    *InSitu_Sz_3D_Seg_Strs_Z_with_Y
    15.0e6,10.0e6,15.0e6
    ~~~   
+ ***<font color=green>InSitu_Sz_3D_Seg_Loca_Z_with_Y</font>** - Locations of Sz stresses when Key_Nonuniform_InSitu_Z_with_Y=1 and Key_InSitu_Strategy = 4.
    Keyword example:
    ~~~bash
    % In this example, along the y-direction, the stress-z in the zone (0.0, 110.0) is 15.0 MPa, the stress-y in the zone (110.0, 140.0) is 10.0 MPa, and the stress-y in the zone (140.0, 250.0) is 15.0 MPa.
    *Key_InSitu_Strategy
    4 
    *Key_Nonuniform_InSitu_Z_with_Y
    1
    *InSitu_Sz_3D_Seg_Strs_Z_with_Y
    15.0e6,10.0e6,15.0e6
    *InSitu_Sz_3D_Seg_Loca_Z_with_Y
    0.0,110.0,140.0,250.0   
    ~~~   
+ ***<font color=green>Key_Gravity</font>** - Gravity.

    = 0, no gravity (default); 

    = 1, apply gravity to the whole model.

    Keyword example:
    ~~~bash
    *Key_Gravity
    1     % Active gravity.      
    ~~~   
+ ***<font color=green>g_X_Y_Z</font>** - Values of gravitational acceleration in $x$, $y$ and $z$ directions.

    1, gravitational acceleration in $x$ direction;

    2, gravitational acceleration in $y$ direction;

    3, gravitational acceleration in $z$ directions.

    Keyword example:
    ~~~bash
    *g_X_Y_Z
    0,0,9.8      % Set gravitational acceleration in z direction.
    ~~~ 
+ ***<font color=green>KEY_RANDOM_SCHEME</font>** - Random generation.

    = 0, randomly generates, but the generated results remain unchanged;

    = 1, randomly generates and the generated results are different (default);

    = 2, randomly generates with a given seed.

    = 3, randomly generates with a seed, seed is randomly generated.

    Keyword example:
    ~~~bash
    *Key_Random_Scheme
    2     
    ~~~ 
+ ***<font color=green>Seed</font>** - Seed for random generation (when *Key_Random_Scheme = 2).

    Keyword example:
    ~~~bash
    *Key_Random_Scheme
    2     
    *Seed
    12345  % Can be arbitrary integer.
    ~~~     
+ ***<font color=green>Key_Crack_Aperture_Method</font>** - Method to calculate crack aperture.

    = 1, equation (default);

    = 2, according to displacements of points on both sides of the fracture surface.

    Keyword example:
    ~~~bash
    *Key_Crack_Aperture_Method
    1     
    ~~~ 

<a id="section-Linear-solvers"></a>

# Solver of linear system of equations 

+ ***<font color=green>Key_SLOE</font>** - Select the solver of linear system.

    = 1, the direct solver (Finding the inverse matrix); 

    = 2, Gaussian elimination;

    = 3, Pardiso (Only available for Intel Fortran compiler); 

    = 4, ITPACK (Iterative method);

    = 5, LAPACK (Default); 

    = 6, MUMPS (Direct method);

    = 7, UMFPACK (Direct method - Multifrontal LU factorization);

    = 8, Lis (Iterative method);

    = 9, SuperLU (Direct method);

    = 10, Y12m (Gaussian elimination);

    = 11, EBE-PCG (Element-by-element preconditioned conjugate gradient solver);

    = 12, STRUMPACK (Direct solver for sparse and dense rank-structured linear systems).

    Keyword example:
    ~~~bash
    *Key_SLOE
    11      % EBE-PCG solver.
    ~~~ 

+ ***<font color=green>Key_K_Sparse</font>** - Assemble the stiffness matrix K in Sparse format (Compressed Sparse Row (CSR) format) to reduce the memory usage. This keyword is available for solvers 6, 7, 8, 9, and 12. If you're working with large-scale models, we recommend either selecting solver 11 or activating sparse matrix storage by setting *Key_K_Sparse = 1 when using other available solvers.
    
    = 0, assemble the stiffness matrix K in full format (default);

    = 1, assemble the stiffness matrix K in Compressed Sparse Row (CSR) format.

    Keyword example:
    ~~~bash
    *Key_SLOE
    8      % Lis solver.
    *Key_K_Sparse
    1      % Assemble the K matrix in CSR format.
    ~~~ 

+ ***<font color=green>Key_EBE_Precondition</font>** - Preconditioner for EBE-PCG solver (solver 11).
    
    = 0, no preconditoner.

    = 1, diagonalized preprocessor (default); 

    = 2, HW preprocessor.

    Keyword example:
    ~~~bash
    *Key_EBE_Precondition
    1    
    ~~~ 

<a id="section-Define-initial-cracks-voids-holes-and-inclusions"></a>

# Define initial cracks, voids (holes) and inclusions

+ ***<font color=green>num_Crack</font>** - Number of initial cracks.

    Keyword example:
    ~~~bash
    *num_Crack
    5     
    ~~~ 
+ ***<font color=green>num_Circ_Hole</font>** - Number of initial circle holes.

    Keyword example:
    ~~~bash
    *num_Circ_Hole
    1
    ~~~     
+ ***<font color=green>num_Circ_Incl</font>** - Number of circular inclusions.

    Keyword example:
    ~~~bash
    *num_Circ_Incl
    5     
    ~~~ 
+ ***<font color=green>num_Poly_Incl</font>** - Number of polygonal inclusions.

    Keyword example:
    ~~~bash
    *num_Poly_Incl
    3     
    ~~~ 
+ ***<font color=green>CRACK_1 - CRACK_50</font>** - Define the coordinates of each initial cracks (line segments) for 2D problem (input format: $P_{1x}$, $P_{1y}$, $P_{2x}$, $P_{2y}$, $P_{3x}$, $P_{3y}$, ... , $P_{nx}$, $P_{ny}$).

    Keyword example:
    ~~~bash
    % Coordinates of 2 initial cracks, the first crack is defined by 2 points, the second crack is defined by 3 points.
    *CRACK_1
    0.0217,0.0229,0.0253,0.0247
    *CRACK_2
    0.0083,0.0193,0.0123,0.0283,0.0223,0.0383
    ~~~ 
+ ***<font color=green>Crack3D_Coor_1 - Crack3D_Coor_50</font>** - Define the coordinates of each initial cracks (plane in 3D composed of four points) for 3D problem (input format: $x_1$, $y_1$, $z_1$, $x_2$, $y_2$, $z_2$, $x_3$, $y_3$, $z_3$, $x_4$, $y_4$, $z_4$).

    Keyword example:
    ~~~bash
    *Crack3D_Coor_1
    10.0,14.0,6.0,10.0,14.0,14.0,10.0,6.0,14.0,10.0,6.0,6.0	
    *Crack3D_Coor_2
    14.0,10.0,6.0,14.0,10.0,14.0,6.00,10.0,14.0,6.0,10.0,6.0
    *Crack3D_Coor_3
    6.0,6.0,10.0,14.0,6.0,10.0,14.0,14.0,10.0,6.0,14.0,10.0
    ~~~ 
+ ***<font color=green>Crack3D_Cir_Coor_1 - Crack3D_Cir_Coor_30</font>** - Define the coordinates of each initial circular cracks for 3D problem (input format: $x$ coordinate of crack center, $y$ coordinate of crack center, $z$ coordinate of crack center, $x$ component of the normal vector of the crack surface, $y$ component of the normal vector of the crack surface, $z$ component of the normal vector of the crack surface, radius $r$ of the crack).

    Keyword example:
    ~~~bash
    *Crack3D_Cir_Coor_1
    8.5,8.5,6.5,1.0,1.0,1.0,2.5
    *Crack3D_Cir_Coor_2
    8.5,8.5,10.5,1.0,1.0,1.0,2.5
    ~~~ 
+ ***<font color=green>Hole_Coor_1 - Hole_Coor_30</font>** - Define the coordinates of initial holes in 2D (input format: $x$, $y$, $r$). 

    Keyword example:
    ~~~bash
    % Define initial holes.
    *Hole_Coor_1
    15.25e-3,7.75e-3,2.0e-3
    *Hole_Coor_2
    4.25e-3,7.75e-3,1.5e-3
    *Hole_Coor_3
    26.25e-3,7.75e-3,1.5e-3
    ~~~ 
+ ***<font color=green>Circ_Inclu_Coor_1 - Circ_Inclu_Coor_30</font>** - Define the coordinates of initial circular inclusions in 2D (input format: $x$, $y$, $r$). 

    Keyword example:
    ~~~bash
    % Define an initial inclusion.
    *Circ_Inclu_Coor_1
    13.5e-3,10.5e-3,3.0e-3
    ~~~ 
+ ***<font color=green>Key_Random_NaCr</font>** - Randomly generate initial cracks or not.

    = 0, do not generate (default);

    = 1, randomly generates;

    = 2, read from FracMan file (*.fab);

    = 3, user-defined cracks (by keywords *Na_CRACK3D_COOR or *NA_CRACK3D_CIR_COOR).

    Keyword example:
    ~~~bash
    % Read initial cracks from FracMan file (*.fab).
    *Key_Random_NaCr
    2
    ~~~     
+ ***<font color=green>Key_NaCr_Active_Scheme_3D</font>** - Activation algorithm for 3D natural fractures.

    = 1, all natural fractures are activated at the initial moment. After being activated by hydraulic fracture, the fractures will be filled with fracturing fluid (default);

    = 2, natural fractures are activated after being communicated by hydraulic fracture. After being activated by hydraulic fracture, the fractures will be filled with fracturing fluid;

    = 3, natural fractures are partially activated after being communicated by hydraulic fracture. After being activated by hydraulic fracture, only part of the crack open and expand along the surface where the natural cracks are located.

    Keyword example:
    ~~~bash
    % Partially enable natural crack when activated by hydraulic cracks.
    *Key_NaCr_Active_Scheme_3D
    3
    ~~~ 
+ ***<font color=green>KIC_NACR</font>** - Fracture toughness $K_{Ic}$ of initial natural cracks.

    Keyword example:
    ~~~bash
    *KIC_NACR
    10
    ~~~     
+ ***<font color=green>KIC_NA_Crack_1 - KIC_NA_Crack_30</font>** - Fracture toughness $K_{Ic}$ of initial natural crack $n$.

    Keyword example:
    ~~~bash
    *KIC_NA_Crack_1
    0.1e6
    *KIC_NA_Crack_2
    0.2e6
    ~~~      
+ ***<font color=green>St_NACR</font>** - Tensile strength $S_t$ of natural cracks.

    Keyword example:
    ~~~bash
    *St_NACR
    0.1e6
    ~~~      
+ ***<font color=green>St_NA_Crack_1 - St_NA_Crack_30</font>** - Tensile strength $S_t$ of natural crack $n$.

    Keyword example:
    ~~~bash
    *St_NA_Crack_1
    0.1e6
    *St_NA_Crack_2
    0.5e6    
    ~~~       
+ ***<font color=green>Key_NaCr_Cross</font>** - Allow the intersection of initially generated natural cracks.
  
    = 1, no (default);

    = 2, yes;

    Keyword example:
    ~~~bash
    *Key_NaCr_Cross
    1
    ~~~                 
+ ***<font color=green>num_Rand_Na_Crack</font>** - Number of cracks need to be randomly generated.

    Keyword example:
    ~~~bash
    *num_Rand_Na_Crack
    10
    ~~~  
+ ***<font color=green>NaCr_Orientation</font>** - Average direction (in degrees) of the cracks need to be randomly generated for 2D problems.

    Keyword example:
    ~~~bash
    *NaCr_Orientation
    -60.0
    ~~~ 
+ ***<font color=green>NaCr_3D_n_Vector</font>** - Average direction (normal vector) of cracks need to be randomly generated for 3D problems.

    Keyword example:
    ~~~bash
    *NaCr_3D_n_Vector 
    1.0,0.0,1.0  
    ~~~     
+ ***<font color=green>NaCr_Ori_Delta</font>** - Fluctuation range (+ or -) of the average direction (in degrees) for 2D problems.

    Keyword example:
    ~~~bash
    *NaCr_Ori_Delta
    10.0
    ~~~ 
+ ***<font color=green>NaCr_3D_n_Vector_Delta</font>** - Amplitude of fluctuation of the normal direction (in degrees) for 3D problems.

    Keyword example:
    ~~~bash
    *NaCr_3D_n_Vector_Delta
    15.0
    ~~~     
+ ***<font color=green>NaCr_Length</font>** - Average length of the cracks need to be randomly generated for 2D problems.

    Keyword example:
    ~~~bash
    *NaCr_Length
    10.0
    ~~~ 
+ ***<font color=green>NaCr_3D_Size</font>** - Average size of the cracks need to be randomly generated for 3D problems.

    Keyword example:
    ~~~bash
    *NaCr_3D_Size   
    10.0
    ~~~     
+ ***<font color=green>NaCr_Len_Delta</font>** - Fluctuation range (+ or -) of the average length of natural cracks for 2D problems.

    Keyword example:
    ~~~bash
    *NaCr_Len_Delta
    2.0
    ~~~ 
+ ***<font color=green>NaCr_3D_Sz_Delta</font>** - Fluctuation range (+ or -) of the average size of natural cracks for 3D problems.

    Keyword example:
    ~~~bash
    *NaCr_3D_Sz_Delta
    5.0
    ~~~  
+ ***<font color=green>num_Rand_Hole</font>** - Number of holes need to be randomly generated.

    Keyword example:
    ~~~bash
    *num_Rand_Hole
    6
    ~~~    
+ ***<font color=green>Rand_Hole_R</font>** - Size of initial generated holes.

    Keyword example:
    ~~~bash
    *Rand_Hole_R
    1.2e-2
    ~~~      
+ ***<font color=green>Rand_Hole_Delta_R</font>** - Fluctuation range (+ or -) of size of initial generated holes.

    Keyword example:
    ~~~bash
    *Rand_Hole_Delta_R
    0.3e-2
    ~~~          
+ ***<font color=green>Key_Hole_Crack_Generate</font>** - Allow cracks to emerge from the inner edge of the holes.

    = 0, no (default);

    = 1, yes. 

    Keyword example:
    ~~~bash
    *Key_Hole_Crack_Generate
    1
    ~~~     
+ ***<font color=green>Key_Num_Cr_Hole_Generated</font>** - The maximum number of cracks emerged from the inner edge of the holes.

    Keyword example:
    ~~~bash
    *Key_Num_Cr_Hole_Generated
    2
    ~~~               
+ ***<font color=green>Key_NaCr_Type_3D</font>** - Shape of initially generated natural cracks for 3D problems.

    = 1, rectangle;

    = 2, circle;

    = 3, polygon.

    Keyword example:
    ~~~bash
    *Key_NaCr_Type_3D
    1
    ~~~         
+ ***<font color=green>Key_Rand_Circ_Incl</font>** - Randomly generate circular inclusions or not.

    = 0, do not generate (default);

    = 1, generates. 

    Keyword example:
    ~~~bash
    *Key_Rand_Circ_Incl
    1
    ~~~ 
+ ***<font color=green>num_Rand_Circ_Incl</font>** - Number of circular inclusions need to be randomly generated.

    Keyword example:
    ~~~bash
    *num_Rand_Circ_Incl
    10
    ~~~ 
+ ***<font color=green>Rand_Circ_Incl_R</font>** - The average radius of the circular inclusions need to be randomly generated.

    Keyword example:
    ~~~bash
    *Rand_Circ_Incl_R
    0.5e-3 
    ~~~ 
+ ***<font color=green>Rand_Circ_Inc_R_Delta</font>** - Fluctuation range (+ or -) of the average radius.

    Keyword example:
    ~~~bash
    *Rand_Circ_Inc_R_Delta
    0.1e-3
    ~~~ 
+ ***<font color=green>Key_Rand_Poly_Incl</font>** - Randomly generate initial regular polygonal inclusions or not.

    = 0, do not generate (default);

    = 1, generates. 

    Keyword example:
    ~~~bash
    *Key_Rand_Poly_Incl
    1
    ~~~ 
+ ***<font color=green>num_Rand_Poly_Incl</font>** - Number of regular polygonal inclusions need to be randomly generated.

    Keyword example:
    ~~~bash
    *num_Rand_Poly_Incl
    10
    ~~~ 
+ ***<font color=green>num_Vert_Poly_Incl</font>** - Number of edges of regular polygonal inclusions.

    Keyword example:
    ~~~bash
    *num_Vert_Poly_Incl
    5
    ~~~ 
+ ***<font color=green>Rand_Poly_Incl_R</font>** - Average radius of circumcircle of regular polygonal inclusions.

    Keyword example:
    ~~~bash
    *Rand_Poly_Incl_R
    0.5e-3 
    ~~~ 
+ ***<font color=green>Rand_Poly_Inc_R_Delta</font>** - Fluctuation range (+ or -) of the average radius of circumcircle of regular polygonal inclusions.

    Keyword example:
    ~~~bash
    *Rand_Poly_Inc_R_Delta
    0.1e-3
    ~~~ 
+ ***<font color=green>Key_User_Defined_2D_Crack_Path</font>** - Enable user defined 2D crack path. The cracks will propagate along user defined paths. This keyword should be used in conjunction with keyword *USER_DEFINED_2D_CRACK_PATH_n. 
    = 1, no (default);

    = 2, yes;

    Keyword example:
    ~~~bash
    *Key_User_Defined_2D_Crack_Path
    1
    ~~~ 
+ ***<font color=green>User_Defined_2D_Crack_Path_1 - User_Defined_2D_Crack_Path_99 </font>** - Define the coordinates of each user defined crack path (Support up to 99 user defined crack paths). Input format of this keyword: $P_{1x}$, $P_{1y}$, $P_{2x}$, $P_{2y}$, $P_{3x}$, $P_{3y}$, ... , $P_{nx}$, $P_{ny}$. This keyword is necessary when keyword *Key_User_Defined_2D_Crack_Path = 1. If multiple crack paths are defined, the program starts simulation from the first path, one by one, until the last crack path.
  
    Keyword example:
    ~~~bash
    % Enable User Defined 2D Crack Path.
    *Key_User_Defined_2D_Crack_Path
    1

    % User Defined 2D Crack Path 1. This crack path consists of 9 coordinate points.
    *User_Defined_2D_Crack_Path_1
    0.154856,0.291413,0.154856,0.480089,0.154856,0.694378,0.253281,0.694577,0.330709,0.670475,0.354331,0.581574,0.333333,0.506059,0.251969,0.484332,0.154856,0.480089

    % User Defined 2D Crack Path 2. This crack path consists of 3 coordinate points.
    *User_Defined_2D_Crack_Path_2
    0.471129,0.290705,0.472441,0.465909,0.473753,0.721977
    ~~~ 

<a id="section-Definition-of-material-parameters"></a>

# Definition of material parameters

+ ***<font color=green>Material_Type_1 - Material_Type_70</font>** - Material type definition for each material.

    = 1, isotropic material;

    = 2, plastic material (Von Mises yield criterion);

    = 3, damage material;

    = 4, Mohr-Coulomb plastic;

    = 5, Composite (for 3D problems only);

    = 6, Drucker-Prager plastic.
 
    Keyword example: 
    ~~~bash
    *MATERIAL_TYPE_1
    1
    *MATERIAL_TYPE_2
    1
    *MATERIAL_TYPE_3
    1
    ~~~
+ ***<font color=green>Material_Para_1 - Material_Para_70</font>** - Definition of material parameters. The following 20 parameters can be defined as needed for each material:

    1, elasticity modulus, $E$;

    2, Poisson's ratio, $ν$;

    3, density, $ρ$;

    4, thickness for 2D plane stress model, $t$;

    5, tensile strength, $σ_t$;

    6, fracture toughness, $K_{Ic}$;

    7, compressive strength, $σ_c$;

    8, coefficient of thermal expansion, $T_α$;

    9, specific heat coefficient, $c$;

    10, conductivity coefficient, $K_{xx}$;

    11, conductivity coefficient, $K_{yy}$;

    12, Mises yield stress $σ_y$;

    13, Tangent modulus after yielding $E_T$;

    14, Hardening parameter $β$;

    15, Cohesive strength $c$;

    16, Internal friction angle $φ$;

    17, Initial damage strain for tension;

    18, Ultimate damage strain for tension;

    19, Initial damage strain for compression;

    20, Ultimate damage strain for compression. 

    Keyword example:
    ~~~bash
    % Parameters for material type 1.
    *Material_Para_1   
    35.0e9,0.25,2000.0,1.0,2.0e6,1.0e6,100.0e6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    % Parameters for material type 2.
    *Material_Para_2   
    0.02e9,0.3,200.0,1.0,1000.0e6,2.0e6,100.0e6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    ~~~
+ ***<font color=green>Material_Para_Added_1 - Material_Para_Added_10</font>** - Additional material parameters required for specific material types (such as material type 5).

    1, elasticity modulus in $x$ direction for material type 5, $E_{1}$;

    2, elasticity modulus in $y$ direction for material type 5, $E_{2}$;

    3, elasticity modulus in $z$ direction for material type 5, $E_{3}$;

    4, Poisson's ratio for material type 5, $ν_{12}$;

    5, Poisson's ratio for material type 5, $ν_{13}$;

    6, Poisson's ratio for material type 5, $ν_{23}$;

    7, Shear elastic modulus for material type 5, $G_{12}$;

    8, Shear elastic modulus for material type 5, $G_{13}$;

    9, Shear elastic modulus for material type 5, $G_{23}$;

    10, Volume fraction of reinforcement phase;

    11, Composite material coordinate system type: 1, Cartesian coordinates; 10, Cylindrical coordinates;

    12, Shear strength, $s_{shear}$;

    13-20, blank.

    Keyword example:
    ~~~bash
    % Example parameters for material type 5.
    *MATERIAL_TYPE_1
    5
    *Material_Para_1   
    38.5e9,0.228,2000.0,1.0,5.0e6,1.0e6,50.0e6,0.5e-5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    *Material_Para_Added_1 
    38.5e9,10.5e9,10.5e9,0.228,0.2,0.2,20e9,5e9,5e9,0.1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    ~~~
+ ***<font color=green>Key_Weibull_E_1 - Key_Weibull_E_10</font>** - Activates Weibull distribution for elastic modulus of material $n$. This keyword must be used in conjunction with WEIBULL_PARAMETERS_E_n.

    = 0, no (default);

    = 1, yes.
 
    Keyword example: 
    ~~~bash
    *Key_Weibull_E_1
    1
    ~~~
+ ***<font color=green>WEIBULL_PARAMETERS_E_1 - WEIBULL_PARAMETERS_E_10</font>** - Weibull parameters of material $n$. 3 parameters are needed to be defined for each material.

    1, shape parameter $k$ (default to 2);

    2, lower bound of the elastic modulus $c0$ (default to 0);

    3, upper bound of the elastic modulus $c1$ (default to infinity). 

    Keyword example:
    ~~~bash
    *Material_Type_1
    1
    *Material_Para_1   
    150.0e9,0.26,7900.0,1.0,1.0e6,1.0e6,100.0e6,0.0,0.0,0.0,0.0,100.0e6,20.0e9,0.5,0.0,0.0,0.0,0.0,0.0,0.0
    *MATERIAL_PARA_ADDED_1
    0,0,0,0,0,0,0,0,0,0,0,300e6,0,0,0,0,0,0,0,0

    % Apply Weibull distribution to the elastic modulus of material 1.
    *Key_Weibull_E_1
    1 
    *Weibull_Parameters_E_1
    2.0,0.0,500e9
    ~~~
<a id="section-Hydraulic-fracturing-simulation"></a>

# Hydraulic fracturing simulation

+ ***<font color=green>Num_Frac</font>** - Number of steps of the hydraulic fracturing simulation.

    Keyword example:
    ~~~bash
    *Num_Frac  
    10
    ~~~
+ ***<font color=green>Key_Symm_HF</font>** - Symmetric hydraulic fracturing model.

    = 0, no (full model);
    
    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Symm_HF
    0
    ~~~
+ ***<font color=green>Cracks_HF_State</font>** - Initial states of each crack (fluid driven or not).

    = 0, free crack without fluid;
    
    = 1, hydraulic fluid-driven crack.

    Keyword example:
    ~~~bash
    *Cracks_HF_State
    1
    ~~~
+ ***<font color=green>Inject_Crack_Num</font>** - For full model (Key_Symm_HF=0), define the crack which contains the injection point of fluid (For symmetric model, the injection point is just the mouth of the edge crack which is also crack number 1). 

    Keyword example:
    ~~~bash
    *Cracks_HF_State
    1
    ~~~
+ ***<font color=green>Inj_Point_Loc</font>** - For full model, define the coordinates of the injection point (input format: $x$, $y$). 

    Keyword example:
    ~~~bash
    *Inj_Point_Loc
    8.5, 8.5
    ~~~
+ ***<font color=green>Inject_Q_Time</font>** - Define the time instants of the data curve of injection rate of fracturing fluid (unit: $s$; 20 time instants at most). 

    Keyword example:
    ~~~bash
    *Inject_Q_Time
    0,10000       
    ~~~
+ ***<font color=green>Inject_Q_Val</font>** - Define the values of injection rate of the data curve of injection rate of fracturing fluid (unit: $m^3/s$).

    Keyword example:
    ~~~bash
    *Inject_Q_Val
    0.001,0.001
    ~~~
+ ***<font color=green>Inject_c_Time</font>** - Define the time instants of the data curve of volumetric concentration of injected proppant (unit: $s$; 20 time instants at most).

    Keyword example:
    ~~~bash
    *Inject_c_Time
    0,10000,20000  
    ~~~
+ ***<font color=green>Inject_c_Val</font>** - Define the values of concentration of the data curve of volumetric concentration of injected proppant.

    Keyword example:
    ~~~bash
    *Inject_c_Val
    0,0.1,0.3
    ~~~
+ ***<font color=green>Key_Visco_Type</font>** - Static viscosity or dynamic viscosity.

    = 1, static viscosity;

    = 2, dynamic viscosity.

    Keyword example:
    ~~~bash
    *Key_Visco_Type
    1
    ~~~
+ ***<font color=green>Viscosity</font>** - The viscosity of fracturing fluid (unit: $Pa·s$).

    Keyword example:
    ~~~bash
    *Viscosity
    0.01
    ~~~
+ ***<font color=green>Viscosity_Par_m</font>** - Parameter $m$ of dynamic viscosity (define when *Key_Visco_Type =2).

    Keyword example:
    ~~~bash
    *Viscosity_Par_m
    2
    ~~~
+ ***<font color=green>Key_Proppant</font>** - Consider proppant or not.

    = 0, no (default);

    = 1, yes. 

    Keyword example:
    ~~~bash
    *Key_Proppant
    1
    ~~~
+ ***<font color=green>Key_Propp_Trans</font>** - Consider the transport of proppant or not.

    = 0, no (default);

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Propp_Trans
    1
    ~~~
+ ***<font color=green>Key_Leakoff</font>** - Consider the leak off of the fracturing fluid or not.

    = 0, no (default);

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Leakoff
    0
    ~~~    
+ ***<font color=green>Coeff_Leak</font>** - Leak coefficient of the Carter model (define when *Key_Leakoff =1).

    Keyword example:
    ~~~bash
    *Coeff_Leak
    1.0d-5
    ~~~    
+ ***<font color=green>Key_Crack_Inner_Pressure</font>** - Active initial fluid pressure of cracks.

    = 0, no (default);

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Crack_Inner_Pressure
    1
    ~~~    
+ ***<font color=green>INI_CRACK_PRESSURE</font>** - Initial fluid pressure of all cracks.

    Keyword example:
    ~~~bash
    *INI_CRACK_PRESSURE
    10.0e6
    ~~~  
+ ***<font color=green>INI_CRACK_PRESSURE_1 - INI_CRACK_PRESSURE_50</font>** - Initial fluid pressure of cracks.

    Keyword example:
    ~~~bash
    *INI_CRACK_PRESSURE_1
    10.0e6
    *INI_CRACK_PRESSURE_2
    10.0e6
    ~~~  
+ ***<font color=green>Key_3D_HF_Time_Step_Method</font>** - Method to solve the fluid-solid coupling problem for the 3D hydraulic fracturing simulation with slick water.

    = 1, Dual-layer Newton-Raphson iteration (default);

    = 2, Bisection.

    Keyword example:
    ~~~bash
    *Key_3D_HF_Time_Step_Method
    1
    ~~~      
+ ***<font color=green>Key_3D_HF_SlipWater_fk_Type</font>** - Calculation method of the equivalent stress intensity factor crack propagation criterion for 3D hydraulic fracturing simulation with slick water.

    = 1, the maximum equivalent stress intensity factor;

    = 2, the avarage equivalent stress intensity factor;

    = 3, the minimum equivalent stress intensity factor.

    Keyword example:
    ~~~bash
    *Key_3D_HF_SlipWater_fk_Type
    3
    ~~~        
+ ***<font color=green>Key_HF_Multistage_3D</font>** -  Active staged fracturing.

    = 0, no (default);

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_HF_Multistage_3D
    1
    ~~~       
+ ***<font color=green>num_Wellbore</font>** -  Number of wellbores.
  
    Keyword example:
    ~~~bash
    *num_Wellbore            
    1 
    ~~~   
+ ***<font color=green>num_Points_WB_1 - num_Points_WB_10</font>** -  Number of points to define wellbore $n$.
  
    Keyword example:
    ~~~bash
    *num_Points_WB_1            
    2 
    ~~~   
+ ***<font color=green>Wellbore_Coors_1_1 - Wellbore_Coors_10_6</font>** -  Points to define wellbore $n$.
  
    Keyword example:
    ~~~bash
    % Point 1 to define wellbore 1.
    *Wellbore_Coors_1_1
    100.0,100.0,0.1
    % Point 2 to define wellbore 1.
    *Wellbore_Coors_1_2
    100.0,100.0,199.9
    ~~~   
+ ***<font color=green>WELLBORES_START_POINT_1 - WELLBORES_START_POINT_10</font>** - Start point for fracturing of wellbore $n$.
  
    Keyword example:
    ~~~bash
    *WELLBORES_START_POINT_1
    100.0,100.0,90.0
    ~~~   
+ ***<font color=green>WELLBORES_END_POINT_1 - WELLBORES_END_POINT_10</font>** - End point for fracturing of wellbore $n$.
  
    Keyword example:
    ~~~bash
    *WELLBORES_END_POINT_1
    100.0,100.0,110.0
    ~~~   
+ ***<font color=green>num_Stages_Wellbores_1 - num_Stages_Wellbores_10</font>** - Number of stages of wellbore $n$.
  
    Keyword example:
    ~~~bash
    *num_Stages_Wellbores_1
    1    
    ~~~   

+ ***<font color=green>NUM_CRS_STAGES_WELLBORES_1_1 - NUM_CRS_STAGES_WELLBORES_5_5</font>** - Number of clusters of a stage of wellbore $n$.
  
    Keyword example:
    ~~~bash
    % Number of clusters of stage 1 of wellbore 5.
    *NUM_CRS_STAGES_WELLBORES_5_1              
    1
    ~~~   
+ ***<font color=green>Key_Gen_Ini_Crack_Wellbores</font>** - Shape of initial hydraulic fractures.

    = 1, rectangle;

    = 2, circle (default);

    = 3, polygon.

    Keyword example:
    ~~~bash
    % Generate circle initial fractures along the wellbore.
    *Key_Gen_Ini_Crack_Wellbores
    2
    ~~~   

+ ***<font color=green>Num_Poly_Edges_PolyCr_WB</font>** - The number of sides (default to 5) of the polygon initial hydraulic fractures when Key_Gen_Ini_Crack_Wellbores = 3.
  
    Keyword example:
    ~~~bash
    % The number of sides (default to 5) of the polygon initial hydraulic fractures.
    *Num_Poly_Edges_PolyCr_WB              
    5
    ~~~ 

+ ***<font color=green>Key_Gen_Ini_Crack_Rule_Wellbores</font>** - Direction of the generated initial hydraulic fractures.

    = 1, the generated initial fractures are perpendicular to the wellbore (default);

    = 2, the generated initial fractures are perpendicular to the direction of the minimum principal stress.

    = 3, user defined direction by giving the normal vector of the generated initial hydraulic fractures (see keyword *Normal_Vector_Gen_Ini_Crack_Wellbores).

    Keyword example:
    ~~~bash
    % Generate circle initial fractures along the wellbore.
    *Key_Gen_Ini_Crack_Wellbores
    2
    % The generated initial fractures are perpendicular to the direction of the minimum principal stress.
    *Key_Gen_Ini_Crack_Rule_Wellbores
    3
    ~~~

+ ***<font color=green>Normal_Vector_Gen_Ini_Crack_Wellbores</font>** - User defined normal vector of the generated initial hydraulic fractures when *Key_Gen_Ini_Crack_Rule_Wellbores=3. The normal vector of fractures need to be defined as: $n_x$, $n_y$, $n_z$.

    Keyword example:
    ~~~bash
    % Generate circle initial fractures along the wellbore.
    *Key_Gen_Ini_Crack_Wellbores
    2
    % User defined direction of the generated initial fractures. The nomal vector of initial fractures in this example is (1,1,1).
    *Normal_Vector_Gen_Ini_Crack_Wellbores
    1,1,1
    ~~~

+ ***<font color=green>Size_Ini_Crack_Wellbores</font>** - Size of initial hydraulic fractures, for circle or polygon cracks, size denotes diameter.

    Keyword example:
    ~~~bash
    *Size_Ini_Crack_Wellbores
    3.0
    ~~~   
+ ***<font color=green>INJECTION_Q_STAGES_WELLBORES_1_1 - INJECTION_Q_STAGES_WELLBORES_5_5</font>** - Injection rate of a stage of wellbore 1 (unit: $m^3/s$).
  
    Keyword example:
    ~~~bash
    % Injection rate of stage 1 of wellbore 5.
    *INJECTION_Q_STAGES_WELLBORES_5_1
    0.01
    ~~~   
+ ***<font color=green>INJECTION_T_STAGES_WELLBORES_1_1 - INJECTION_T_STAGES_WELLBORES_5_5</font>** - Injection time of a stage of wellbore 1 (unit: $s$).
  
    Keyword example:
    ~~~bash
    % Injection time of stage 1 of wellbore 5.
    *INJECTION_T_STAGES_WELLBORES_5_1
    1000
    ~~~   
+ ***<font color=green>Key_Get_Permeability</font>** - Calculate the element equivalent permeability for 3D hydraulic fracturing simulation.

    = 0, no;

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_Get_Permeability
    1
    ~~~  
+ ***<font color=green>Key_3D_Slip_HF_Keep_Pressure</font>** - Keep the fluid pressure after each stage for 3D slip water hydraulic fracturing simulation.

    = 0, no;

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_3D_Slip_HF_Keep_Pressure
    1
    ~~~
    

<a id="section-Nonlinear-analysis"></a>

# Nonlinear analysis

+ ***<font color=green>NL_ITRA</font>** - Maximum number of Newton-Raphson iteration (default：30).
  
    Keyword example:
    ~~~bash
    *NL_ITRA
    40
    ~~~    
+ ***<font color=green>NL_ATOL</font>** - Maximum Norm-2 value of the residual of Newton-Raphson iteration (default：1.0e8).

    Keyword example:
    ~~~bash
    *NL_ATOL
    2.0e8
    ~~~    
+ ***<font color=green>NL_NTOL</font>** - Maximum number of bisection of force for Newton-Raphson iteration (default：6).

    Keyword example:
    ~~~bash
    *NL_NTOL
    10
    ~~~    
+ ***<font color=green>NL_TOL</font>** - Convergence tolerance for Newton-Raphson iteration (default：1.0e-6).

    Keyword example:
    ~~~bash
    *NL_TOL
    1.0e6
    ~~~    
+ ***<font color=green>NL_TIMS_1 - NL_TIMS_30</font>** - Load step control, 5 parameters are needed for each load step:

    1, starting time for the current load step;

    2, ending time for the current load step;

    3, time increment for the current load step;

    4, starting force factor for the current load step;

    5, ending force factor for the current load step.

    Keyword example:
    ~~~bash
    *NL_TIMS_1
    0.0,0.4,0.1,0.0,0.5
    *NL_TIMS_2
    0.4,1.0,0.05,0.5,1.0
    ~~~    

<a id="section-Cohesive-crack"></a>

# Cohesive crack

+ ***<font color=green>Coh_Constitutive_type</font>** - Constitutive model of the cohesive crack: 

    = 1, bilinear model, first rise and then fall (Define both *Coh_Width_Critical1 and *Coh_Width_Critical2); 

    = 2, linear model (Define *Coh_Width_Critical2); 

    = 3, constant model (Define *Coh_Width_Critical2). 

    Keyword example:
    ~~~bash
    *Coh_Constitutive_type
    1
    ~~~    
+ ***<font color=green>Coh_Width_Critical1</font>** - Normal width of crack at which the crack surface has the ultimate normal traction (needs to be defined only when *Coh_Constitutive_type =1).

    Keyword example:
    ~~~bash
    *Coh_Width_Critical1
    2.0e-4
    ~~~    
+ ***<font color=green>Coh_Width_Critical2</font>** - Normal width of crack at which the crack surface has no normal traction.

    Keyword example:
    ~~~bash
    *Coh_Width_Critical2
    1.0e-5
    ~~~    
+ ***<font color=green>Coh_f_Ultimate</font>** - The ultimate normal traction.

    Keyword example:
    ~~~bash
    *Coh_f_Ultimate
    3.0D5
    ~~~    
+ ***<font color=green>Coh_Tangential_Key</font>** - Consider tangential traction or not.

    = 0, no (default);

    = 1, yes (Define *Coh_Width_Critical1_T, *Coh_Width_Critical2_T, *Coh_f_Ultimate_T).

    Keyword example:
    ~~~bash
    *Coh_Tangential_Key
    1
    ~~~    
+ ***<font color=green>Coh_Width_Critical1_T</font>** - Tangential width of crack at which the crack surface has the ultimate tangential traction (needs to be defined only when *Coh_Constitutive_type =1 and * Coh_Tangential_Key =1).

    Keyword example:
    ~~~bash
    *Coh_Width_Critical1_T
    0.25D-3
    ~~~    
+ ***<font color=green>Coh_Width_Critical2_T</font>** - Tangential width of crack at which the crack surface has no tangential traction.

    Keyword example:
    ~~~bash
    *Coh_Width_Critical2_T
    0.50D-3
    ~~~    
+ ***<font color=green>Coh_f_Ultimate_T</font>** - The ultimate tangential traction

    Keyword example:
    ~~~bash
    *Coh_f_Ultimate_T
    5.0D6
    ~~~    

<a id="section-Dynamic-analysis"></a>

# Dynamic analysis

+ ***<font color=green>IDy_Num_Iteras</font>** - Number of steps of the implicit dynamic analysis.

    Keyword example:
    ~~~bash
    *IDy_Num_Iteras
    200
    ~~~    
+ ***<font color=green>IDy_Num_force_Itr</font>** - Number of steps with force applied (In other words, the applied force will be removed if the current step number is larger than *IDy_Num_force_Itr) for the implicit dynamic analysis.

    Keyword example:
    ~~~bash
    *IDy_Num_force_Itr
    100
    ~~~    
+ ***<font color=green>Delt_Time_NewMark</font>** - Time increment of the Newmark time integration algorithm (default to 1.0e-6) for the implicit dynamic analysis. 

    Keyword example:
    ~~~bash
    *Delt_Time_NewMark
    5.0e-6
    ~~~    
+ ***<font color=green>edy_num_iteras</font>** - Number of steps of the explicit dynamic analysis.

    Keyword example:
    ~~~bash
    *edy_num_iteras
    200
    ~~~    
+ ***<font color=green>EDy_num_force_itr</font>** - Number of steps with force applied (In other words, the applied force will be removed if the current step number is larger than *edy_num_iteras) for the explicit dynamic analysis.

    Keyword example:
    ~~~bash
    *EDy_num_force_itr
    100
    ~~~    
+ ***<font color=green>delt_time_explicit</font>** - Time increment of the time step (default to 1.0e-6) for the explicit dynamic analysis. 

    Keyword example:
    ~~~bash
    *delt_time_explicit
    5.0e-6
    ~~~        
+ ***<font color=green>Factor_Prop_Dy</font>** - Factor of dynamic propagation length of cracks (default to 1.63). Propagation length $\Delta l=*Factor_Prop_Dy \times l_c$, where $l_c$ represents the average size of enriched elements. 

    Keyword example:
    ~~~bash
    *Factor_Prop_Dy
    1.63
    ~~~        
+ ***<font color=green>Key_EQ</font>**: Earthquake analysis.

    = 0, no (default); 

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_EQ
    1
    ~~~    
+ ***<font color=green>num_EQ_Ac_nodes</font>** - Number of the nodes with earthquake acceleration applied.

    Keyword example:
    ~~~bash
    *num_EQ_Ac_nodes
    31
    ~~~    
+ ***<font color=green>EQ_Ac_nodes</font>** - List of nodes with earthquake acceleration applied.

    Keyword example:
    ~~~bash
    *EQ_Ac_nodes
    187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217
    ~~~    
+ ***<font color=green>EQ_Ac_Time_Gap</font>** - Time interval of the earthquake acceleration data (define when *Key_EQ=1).

    Keyword example:
    ~~~bash
    *EQ_Ac_Time_Gap
    0.01
    ~~~    

<a id="section-Coupling-of-degrees-of-freedom"></a>

# Coupling of degrees of freedom

+ ***<font color=green>num_CP_x_nodes</font>** - The total amount of nodes need to be coupled in $x$ direction.

    Keyword example:
    ~~~bash
    *num_CP_x_nodes
    21
    ~~~    
+ ***<font color=green>CP_x_nodes</font>** - The nodes number of all the nodes need to be coupled in $x$ direction.

    Keyword example:
    ~~~bash
    *CP_x_nodes
    2,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101
    ~~~    
+ ***<font color=green>num_CP_y_nodes</font>** - The total amount of nodes need to be coupled in $y$ direction.

    Keyword example:
    ~~~bash
    *num_CP_y_nodes
    21
    ~~~    
+ ***<font color=green>CP_y_nodes</font>** - The nodes number of all the nodes need to be coupled in $y$ direction.

    Keyword example:
    ~~~bash
    *CP_y_nodes
    2,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101
    ~~~    

<a id="section-Deactivation-elements"></a>

# Deactivation elements

+ ***<font color=green>Key_EKILL</font>** - Deactivate elements during simulation.

    = 1, no (default);

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_EKILL
    1
    ~~~  
+ ***<font color=green>Ele_Killed_Each_Load_Step_02 - Ele_Killed_Each_Load_Step_40</font>** - Define element lists for load step $n$.

    Keyword example:
    ~~~bash
    % Define elements to be killed at each step (from step 2 to step 5).
    *Ele_Killed_Each_Load_Step_02
    501,466
    *Ele_Killed_Each_Load_Step_03
    502,467  
    *Ele_Killed_Each_Load_Step_04 
    503,468   
    *Ele_Killed_Each_Load_Step_05
    504,469  
    ~~~          

<a id="section-Thermal-load-prestress"></a>

# Thermal load and prestress

+ ***<font color=green>KEY_THERMAL_STRESS</font>** - Enable thermal stress.

    = 0, no (default);

    = 1, yes.

    Keyword example:
    ~~~bash
    *KEY_THERMAL_STRESS
    1
    ~~~  

+ ***<font color=green>Key_Initial_Temperature</font>** - Methods to define initial temperature.

    = 1, define initial temperature according to material type number (default);

    = 2, read initial temperature from file.

    Keyword example:
    ~~~bash
    *Key_Initial_Temperature
    1
    ~~~  

+ ***<font color=green>THERMAL_STR_TEMPER_1 - THERMAL_STR_TEMPER_25</font>** - Initial temperature of material type $n$.

    Keyword example:
    ~~~bash
    *THERMAL_STR_TEMPER_1
    70.0
    *THERMAL_STR_TEMPER_2
    100.0    
    ~~~ 

+ ***<font color=green>Key_Scheme_Thermal_Stress</font>** - Methods to calculate thermal stress.

    = 1, thermal stress is calculated directly from the temperature value (default);

    = 2, thermal stress is calculated based on temperature difference.

    Keyword example:
    ~~~bash
    *Key_Scheme_Thermal_Stress
    1
    ~~~    

+ ***<font color=green>Key_InStress_for_Mat</font>** - Apply prestress to a material.

    = 0, no (default);

    = 1, yes.

    Keyword example:
    ~~~bash
    *Key_InStress_for_Mat
    1
    ~~~  

+ ***<font color=green>Mat_Number_of_InStress_1 - Mat_Number_of_InStress_10</font>** - Material type number to apply prestress.

    Keyword example:
    ~~~bash
    *Mat_Number_of_InStress_1 
    2     % Apply prestress to material type 2.   
    ~~~ 

+ ***<font color=green>MAT_INSTRESS_X, MAT_INSTRESS_Y, MAT_INSTRESS_Z</font>** - Prestress in $x$, $y$, and $z$ directions for a specified material type.

    Keyword example:
    ~~~bash
    *MAT_INSTRESS_X  
    10.0e6
    *MAT_INSTRESS_Y  
    5.0e6
    *MAT_INSTRESS_Z  
    5.0e6
    ~~~ 

<a id="section-Surface-load"></a>

# Surface load

+ ***<font color=green>Num_Surface_Loads</font>** - Number of applied surface load.

    Keyword example:
    ~~~bash
    *Num_Surface_Loads
    3
    ~~~  
+ ***<font color=green>File_Surface_Load_1 - File_Surface_Load_20</font>** -  File suffix of surface load $n$.

    Keyword example:
    ~~~bash
    *File_Surface_Load_1
    surf
    *File_Surface_Load_2
    top    
    *File_Surface_Load_3
    inner       
    % Each line of the *.surf file stores surface element number, node 1, node 2, node 3, node 4, and the area of the surface element, for example:
    % 1415.    7789.    7828.    7790.    7829.   0.12500000E-06
    % 1416.    7828.    7867.    7829.    7868.   0.12500000E-06
    ~~~ 
+ ***<font color=green>Surface_Pressure_1 - Surface_Pressure_20</font>** - Pressure value if surface load $n$.

    Keyword example:
    ~~~bash
    *Surface_Pressure_1
    -10.0e6
    ~~~ 

<a id="section-Control-of-the-program"></a>

# Control of the program

+ ***<font color=green>Key_Clear_All</font>** - Clear all the results files in the current work directory before starting.

    = 0, disable (default);

    = 1, enable.

    Keyword example:
    ~~~bash
    *Key_Clear_All
    1
    ~~~    
+ ***<font color=green>Key_Close_Window</font>** - Close the console window when calculation is finished.

    = 0, wait user to close (default); 

    = 1, auto close.

    Keyword example:
    ~~~bash
    *Key_Close_Window
    0
    ~~~    
+ ***<font color=green>Key_Data_Format</font>** - Save data in ASCII format or in binary format.

    = 1, ASCII format (default); 

    = 2, binary format.

    Keyword example:
    ~~~bash
    *Key_Data_Format
    1
    ~~~    
+ ***<font color=green>Key_Num_Process</font>** - Set the number of threads of CPU for OpenMP (if taken as 99, then all threads are available, default to 1).

    Keyword example:
    ~~~bash
    *Key_Num_Process
    4     % Use 4 threads.
    ~~~   
+ ***<font color=green>Key_Save_vtk</font>** - Save the simulation results as vtk files for postprocessing using ParaView.

    = 0, do not save vtk file; 

    = 1, save vtk file (default).

    Keyword example:
    ~~~bash
    *Key_Save_vtk
    1
    ~~~      
+ ***<font color=green>Key_Simple_Post</font>** - Save only necessary result files.

    = 0, save all result files (default); 

    = 1, save only necessary result files.

    Keyword example:
    ~~~bash
    *Key_Simple_Post
    1
    ~~~       
+ ***<font color=green>Key_Post_CS_N_Strs</font>** - Calcutale and save stress results of nodes.

    = 0, no; 

    = 1, yes (default).

    Keyword example:
    ~~~bash
    *Key_Post_CS_N_Strs
    0
    ~~~          
+ ***<font color=green>Key_Play_Sounds</font>** - Play sound using Python script when finised the simulation.

    = 0, no (default); 

    = 1, yes.

    Keyword example:
    ~~~bash
    % Play sound (using python).
    *Key_Play_Sounds
    1
    ~~~  

<a id="section-Keywords-file-examples"></a>

# Keywords file examples    

+ 2D FEM analysis.
    ~~~bash
    % This is a simple PhiPsi kpp file example.
    % Pay attention to *Work_Directory, some modification might be necessary!

    % Working directory.
    *Work_Directory
    X:\PhiPsi Work\FEM

    % Filename of input files.
    *Filename
    FEM
        
    % Parameters.
    Elasticity_modulus = 70.0e9
    Poisson_ratio      = 0.3

    % Analysis type (Quasi-static).
    *Key_Analysis_Type
    1

    % Plane stress.
    *Key_Type_2D
    1

    % Linear system solver (SuperLU).
    *Key_SLOE
    9

    % Number of propagation steps.
    *Num_Substeps
    1

    % Material(1-E,2-v,3-density,4-thick,5-St,6-KIc,7-Sc,8-20(blank))
    *Material_Para_1
    Elasticity_modulus,Poisson_ratio,2700.0,1.0,1.0e6,1.0e6,100.0e6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

    % Play sound (using python).
    *Key_Play_Sounds
    1
    ~~~ 

+ 2D XFEM analysis.

    ~~~bash
    % Working directory.
    *Work_Directory
    X:\PhiPsi Work\XFEM

    % Filename of input files.
    *Filename
    XFEM

    % Number of initial cracks.
    *num_Crack
    1

    % Analysis type (Quasi-static).
    *Key_Analysis_Type
    1

    % Plane strain.
    *Key_Type_2D
    2

    % Calculate stress intensity factors using the interaction integral method.
    *Key_SIFs_Method
    2

    % Linear system solver.
    *Key_SLOE
    % 5 %LAPACK
    11 %PCG-EBE

    % Number of propagation steps.
    *Num_Substeps
    1

    % Material(1-E,2-v,3-density,4-thick,5-St,6-KIc,7-Sc,8-20(blank))
    *Material_Para_1   
    20.0e9,0.2,2000.0,1.0,1.0e6,1.0e6,100.0e6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

    % Coordinates of initial cracks.
    *CRACK_1
    0.00615,0.013,0.01485,0.018

    % Play sound (using python).
    *Key_Play_Sounds
    1
    ~~~ 


+ 3D XFEM analysis.

    ~~~bash
    % Working directory.
    *Work_Directory
    X:\PhiPsi_Work\exa_3D_block_tension_cir

    % Finename of input files.
    *Filename
    exa_3D_block_tension_cir

    % Number of CPU threads.
    *Key_Num_Process
    4

    % 3-D problem.
    *Key_Dimension
    3

    % Number of initial cracks.
    *num_Crack
    2

    % Analysis type (Quasi-static).
    *Key_Analysis_Type
    1

    % Linear system solver.
    *Key_SLOE
    11     % PCG-EBE.

    % Number of propagation step.
    *Num_Substeps
    5

    % Crack propagation criterion.
    *CFCP
    2     % Weighted average maximum principal tensile stress criterion.

    % Material
    *Material_Para_1   
    20.0e9,0.3,2000.0,1.0,0.1e6,1.0e6,100.0e6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

    % Define the coordinates of each initial cracks for 3D problem.
    % (input format: x_of_circle_center, y_of_circle_center, z_of_circle_center, x_of_direction_vector, y_of_direction_vector, z_of_direction_vector, radius)
    *Crack3D_Cir_Coor_1
    8.5,8.5,6.5,1.0,1.0,1.0,2.5
    *Crack3D_Cir_Coor_2
    8.5,8.5,10.5,1.0,1.0,1.0,2.5

    % Allow crack 1 and 2 to propagate.
    *Cracks_Allow_Propa
    1,1

    % Do not close the DOS Window when the calculation is finished (default).
    *Key_Close_Window
    0

    % Save crack radius
    *Key_Save_Crack_Radius
    1
    ~~~ 