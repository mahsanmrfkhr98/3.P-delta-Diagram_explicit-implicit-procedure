# Clear Previous Data
wipe
# units: lb, psi, in
# Shekl, safheye 2 az file "Session 03.pdf".
# Define Model Properties: model  BasicBuilder  -ndm  Number_of_Dimensions  -ndf  Nunber_of_DoFs_in_Each_Node
						   model  BasicBuilder  -ndm         2              -ndf              2
						   
# Define Nodes Information
	# Create nodes: node  node#   x[in]  y[in]
					node    1      0.0    0.0
					node    2    160.0    0.0
					node    3    320.0    0.0
	
	# Assign nodes to a Group: region  Region#  -Type  Target_Nodes 
				   			   region    1      -node      1 2 3
	
	# Set Restrained DoFs: fix  Node#  X_Restrained?  Y_Restrained?
						   fix    1         1                1
					       fix    3         1                1
# Define Elements Information
    # Create Material with Initial Stress: uniaxialMaterial  Steel02  Material#  Fy(Anything)[psi]    E[psi]    alpha(Linear)[-]        Some_Default_Parameters          Initial_Stress(N0/A)[psi] 
	#                                                                                                                             -------------------------------------
	#                                                                                                                              R0    cR1   cR2   a1   a2   a3   a4
						                   uniaxialMaterial  Steel02     1            100.0         36000000.0       1.0          20.0  0.925  0.15  0.0  1.0  0.0  1.0            300
	# Create Corotational Truss Elements:    element  corotTruss  Element#  Start_Node  End_Node  A[in^2]  Assigned_Material#
							                 element  corotTruss     1           1         2        5.0            1
							                 element  corotTruss     2           2         3        5.0            1	
	# Assign Elements to a Group:  region  Region#  -Type  Target_Elements 
				   			       region    2      -ele        1 2 
	
# Define Applied Loads
	# Create a Linear Loading TimeSeries: timeSeries  Linear  TimeSeries#
								          timeSeries  Linear      1
									 
	# Create a Plain load pattern associated with the TimeSeries: pattern  Plain  pattern#  Assigned_TimeSeries#  {
																  pattern  Plain     1              1             {
	# Create the nodal loads:																					     load  Node#  X_Value[lb]  Y_Value[lb]
																												     load    2       0.0          0.1
																												  }

# Define Analysis Parameters																							  
							system BandGeneral
							numberer RCM
							constraints Transformation 						
							integrator LoadControl 1.0					
                            test NormUnbalance 0.000000001 20000 0
							algorithm Newton
							analysis Static
							
# Create Output Files
	# Output File for Nodal Info:    recorder  Node     -file   Output_File_Name          -time[?]  -node/Region  Node/Region#  -dof  Target_DoFs  Type
								     recorder  Node     -file   OSNodalDisplacements.txt            -region           1         -dof      1 2      disp   
									 
# Start the Analysis: analyze  Number_of_Steps									 
					  analyze       2000										 