# Clear Previous Data
wipe
# units: lb, psi, in
# Shekl, safheye 2 az file "Session 03.pdf".
# Define Model Properties: model  BasicBuilder  -ndm  Number_of_Dimensions  -ndf  Nunber_of_DoFs_in_Each_Node
						   model  BasicBuilder  -ndm         2              -ndf              3
						   
# Define Nodes Information
	# Create nodes: node  node#   x[in]  y[in]
					node    1      0.0    0.0
					node    2    160.0    0.0
					node    3    320.0    0.0

	# Assign actual nodes to a Group: region  Region#  -Type  Target_Nodes 
				   			          region    1      -node      1 2 3

	# Create virtual nodes with the same coordinates as the actual nodes: node  node#(...1 means right, ..2 means left)  x[in]  y[in]
					                                                      node                  11                         0.0   0.0
					                                                      node                  21                       160.0   0.0
																		  node                  22                       160.0   0.0
					                                                      node                  32                       320.0   0.0
						   
	# Set Restrained DoFs: fix  Node#  X_Restrained?  Y_Restrained?  Rotation_Restrained?
						   fix    1         1               1                1
					       fix    3         1               1                1
	# Restrain the virtual nodes to the actual nodes in X and Y directions: equalDOF  Master_Node#  Slave_Node#  Preferred_DoFs 
                                                                            equalDOF      1            11             1 2
                                                                            equalDOF      2            22             1 2																			
                                                                            equalDOF      2            21             1 2																			
                                                                            equalDOF      3            32             1 2
	# Define Elements Information
    # Create Material with Initial Stress: uniaxialMaterial  Steel02  Material#  Fy(Anything)[psi]    E[psi]    alpha(Linear)[-]        Some_Default_Parameters          Initial_Stress(N0/A)[psi] 
	#                                                                                                                             -------------------------------------
	#                                                                                                                              R0    cR1   cR2   a1   a2   a3   a4
						                   uniaxialMaterial  Steel02     1            100.0         36000000.0       1.0          20.0  0.925  0.15  0.0  1.0  0.0  1.0            300
    # Create Cross-Section for Beam-Column Elements:  section  Fiber  Section#  {  patch  Type(circ/rect/quad)  Assigned_Material#  Number_of_Fibers_in_Perimeter  Number_of_Fibers_in_Radius  Center_Y[in]  Center_Z[in]  Internal_Rdius[in]  External_Rdius[in]  Start_Angle[degree]  End_Angle[degree] }
                                                      section  Fiber     1      {  patch         circ                   1                          2                              2                 0.0          0.0              0.0               1.261                0.0                360.0        }
    # Create Geometric Transformer for Beam-Column Elements:    geomTransf  Type[Linear/PDelta]  Transformer#	
	                                                            geomTransf              PDelta        1
                                                                geomTransf              PDelta        2
	
	# Create Beam-Column Elements:                    element  forceBeamColumn  Element#  Left_Node  Right_Node  Number_of_Integration_Points  Section#  Transformer#
                                                      element  forceBeamColumn      1        11          22                    3                  1           1
                                                      element  forceBeamColumn      2        21          32                    3                  1           2

	# Assign Elements to a Group:  region  Region#  -Type  Target_Elements 
				   			       region    2      -ele        1 2                                                                  										   
	# Create Elastic Material for Zero-Length Springs:  uniaxialMaterial  Elastic  Material#  E(Very_Low)[psi]
						                                uniaxialMaterial  Elastic     11      0.00000000000001
						                                uniaxialMaterial  Elastic     22      0.00000000000001
						                                uniaxialMaterial  Elastic     21      0.00000000000001
						                                uniaxialMaterial  Elastic     32      0.00000000000001
														
	#												                                                                                                                     http://opensees.berkeley.edu/wiki/index.php/ZeroLength_Element
	# Define Zero-Length Spring Elements Between the Actual and the virtual Nodes:  element  zeroLength  Element#  Start_Node  End_Node  -mat  Assigned_Material#  -dir  Preffered_Material_Direction												
                                                                                    element  zeroLength     11          1         11     -mat          11          -dir             6
                                                                                    element  zeroLength     22          2         22     -mat          22          -dir             6
                                                                                    element  zeroLength     21          2         21     -mat          21          -dir             6
                                                                                    element  zeroLength     32          3         32     -mat          32          -dir             6
# Define Applied Loads
	# Create a Linear Loading TimeSeries: timeSeries  Linear  TimeSeries#
								          timeSeries  Linear      1
									 
	# Create a Plain load pattern associated with the TimeSeries: pattern  Plain  pattern#  Assigned_TimeSeries#  {
																  pattern  Plain     1              1             {
	# Create the nodal loads:																					     load  Node#  X_Value[lb]  Y_Value[lb]  Moment_Value[lb]
																												     load    2       0.0          0.1            0.0
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
					  analyze      292																				