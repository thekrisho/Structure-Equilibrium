# Structure-Equilibrium

This experiment was conducted to explore truss structure equilibrium by writting a universal program that could solve any two dimensional (N x N) static truss using the Method of Joints. 

Use:
- Create a text file (see Sample Truss.txt) containing member-joint connectivity, the coordinates of each joint, the reaction forces at each joint and all of the external forces acting on the truss. Our program will prompt you for the name of the text file containing the truss information. 
- The data will then be placed into arrays from which it will construct the matrices M and E, which contain all of the necessary information on the truss structure to compute the external forces. 
- Each member force will then be computed and printed in the screen, as well as written to an output file. 

Background:
- Referring to the truss sample, we can begin analyzing the joints coordinates and their connectivity to the members. We then observe the external forces and the different supports the structure can handle.
- Once this is done we can check to make sure the system is solvable. We have 4 joints each with forces in X and Y directions which gives us 8 equations in total. We have 5 unknown member forces as well as 3 unknown reaction forces. This system is solvable. Now we can jump into solving the sum of forces in X and Y for each and every individual joint. 


                                                          ΣFx=0
                                                          ΣFy=0
                                                          
- Using the coordinates of each joint we can calculate the length of each members as well as theta angles. 

                                            Li = SQRT[ (xjj-xj)^2 + (yjj-yj)^2 ] 
                                            cos (θi) = (xjj-xj/Li) + (xto-xfrom/Li)
                                            sin (θi) = (yjj-yj/Li) + (yto-yfrom/Li)
                                            
- Now that we have all acquired all the data to build this truss in 2D space we can put this system into matrix form.


