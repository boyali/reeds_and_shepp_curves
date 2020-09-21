'''
   ########## ----------------------------- CURVES START HERE ----------------------------------------
   Every curve family has four type of equation. Usually only one equation is given and the others are derived using
   the "Time Flip" or Reflection operations.

   Time flip is simply reversing all the curves. It's derivation is easy. If you flip theta -> - theta, the x and y
   changes according to trigonometry
'''

path_dict = {}  # num: maneuvers

'''
  ########## ---------------------------------- C|C|C Curves --------------------------------------------   
  Curve | Curve | Curve Family  - Curves from [1 -> 4]   "|" cusp means direction change between the curves

  num 1     : ["L+", "R-", "L+"], arguments -- c_c_c(x,y,phi,ap,b1)     -> is given     
  num 2     : ["L-", "R+", "L-"], arguments -- c_c_c(-x,y,-phi,am,b1)   -> time flip of num1 
  num 3     : ["R+", "L-", "R+"], arguments -- c_c_c(x,-y,-phi,am,b1)   -> reflection of num1 
  num 4     : ["R-", "L+", "R-"], arguments -- c_c_c(-x,-y,phi,ap,b1)   -> time flip + reflection of num 1   


'''

path_dict[1] = ["L+", "R-", "L+"]
path_dict[2] = ["L-", "R+", "L-"]
path_dict[3] = ["R+", "L-", "R+"]
path_dict[4] = ["R-", "L+", "R-"]

'''
  ########## ---------------------------------- C|CC Curves --------------------------------------------   


  num 5     : ["L+", "R-", "L-"], arguments -- c_cc(x,y,phi,ap,b1)        -> is given
  num 6     : ["L-", "R+", "L+"], arguments -- c_cc(-x,y,-phi,am,b1)      -> time flip of num 5
  num 7     : ["R+", "L-", "R-"], arguments -- c_cc(x,-y,-phi,am,b1)     -> reflection of num 5
  num 8     : ["R-", "L+", "R+"], arguments -- c_cc(-x,-y,phi,ap,b1)      -> time flip and reflection of num 5 

'''

path_dict[5] = ["L+", "R-", "L-"]
path_dict[6] = ["L-", "R+", "L+"]
path_dict[7] = ["R+", "L-", "R-"]
path_dict[8] = ["R-", "L+", "R+"]

'''
  ########## ---------------------------------- CSC Curves --------------------------------------------   
  Curve - Straight - Curve Family  - Curves from [9 -> 16]

  num 9     : ["L+", "S+", "L+"], arguments -- csca(x,y,phi,ap,b1)      -> given 
  num 10    : ["R+", "S+", "R+"], arguments -- csca(x,-y,-phi,am,b1)    -> reflection of num 9 
  num 11    : ["L-", "S-", "L-"], arguments -- csca(-x,y,-phi, am, b1)  -> time flip of num 9
  num 12    : ["R-", "S-", "R-"], arguments -- csca(-x,-y,phi,ap,b1)    -> time flip + reflection num9

  num 13    : ["L+", "S+", "R+"], arguments -- cscb(x,y,phi,ap,b2)      -> given
  num 14    : ["R+", "S+", "L+"], arguments -- cscb(x,-y,-phi,am,b2)    -> reflection of num 13
  num 15    : ["L-", "S-", "R-"], arguments -- cscb(-x,y,-phi,am,b2)    -> time flip of num 13
  num 16    : ["R-", "S-", "L-"], arguments -- cscb(-x,-y,phi,ap,b2)    -> time flip and reflection of num 13  
'''
path_dict[9] = ["L+", "S+", "L+"]
path_dict[10] = ["R+", "S+", "R+"]
path_dict[11] = ["L-", "S-", "L-"]
path_dict[12] = ["R-", "S-", "R-"]

path_dict[13] = ["L+", "S+", "R+"]
path_dict[14] = ["R+", "S+", "L+"]
path_dict[15] = ["L-", "S-", "R-"]
path_dict[16] = ["R-", "S-", "L-"]
'''
  ########## ---------------------------------- C Cu | Cu C  Curves -------------------------------------------------   
  Curves from [17 -> 20]
  num 17     : ["L+", "R+", "L-", "R-"], arguments -- cu_cuc(x,y,phi,ap,b2)         -> is given
  num 18     : ["R+", "L+", "R-", "L-"], arguments -- ccu_cuc(x,-y,-phi,am,b2)      -> reflection of num 17    
  num 19     : ["L-", "R-", "L+", "R+"], arguments -- ccu_cuc(-x,y,-phi,am,b2)      -> time filip of num 17   
  num 20     : ["R-", "L-", "R+", "L+"], arguments -- ccu_cuc(-x,-y,phi,ap,b2)      -> time flip and reflect of num 17 

'''
path_dict[17] = ["L+", "R+", "L-", "R-"]
path_dict[18] = ["R+", "L+", "R-", "L-"]
path_dict[19] = ["L-", "R-", "L+", "R+"]
path_dict[20] = ["R-", "L-", "R+", "L+"]

'''
 ########## ---------------------------------- C | Cu Cu | C     Curves -------------------------------------------- 
 num 21     : ["L+", "R-", "L-", "R+"], arguments -- c_cucu_c(x,y,phi,ap,b2)   --> given
 num 22     : ["R+", "L-", "R-", "L+"], arguments -- c_cucu_c(x,-y,-phi,am,b2) --> reflection of num 21 
 num 23     : ["L-", "R+", "L+", "R-"], arguments -- c_cucu_c(-x,y,-phi,am,b2) --> time flip of num 21
 num 24     : ["R-", "L+", "R+", "L-"], arguments -- c_cucu_c(-x,-y,phi,ap,b2) --> time flip and reflection

'''
path_dict[21] = ["L+", "R-", "L-", "R+"]
path_dict[22] = ["R+", "L-", "R-", "L+"]
path_dict[23] = ["L-", "R+", "L+", "R-"]
path_dict[24] = ["R-", "L+", "R+", "L-"]

'''
 ########## ---------------------------------- C | C[pi/2] S C    Curves -------------------------------------------- 
 num 25     : ["L+", "R-", "S-", "L-"], arguments --  c_c2sca(x,y,phi,ap,b1)        -> given 
 num 26     : ["R+", "L-", "S-", "R-"], arguments --  c_c2sca(x,-y,-phi,am,b1)      -> reflection of num 25
 num 27     : ["L-", "R+", "S+", "L+"], arguments --  c_c2sca(-x,y,-phi,am,b1)      -> time flip of num 25
 num 28     : ["R-", "L+", "S+", "R+"], arguments --  c_c2sca(-x,-y,phi,ap,b1)      -> time flip of num 26

 num 29     : ["L+", "R-", "S-", "R-"], arguments --  c_c2scb(x,y,phi,ap,b2)        -> given
 num 30     : ["R+", "L-", "S-", "L-"], arguments --  c_c2scb(x,-y,-phi,am,b2)      -> reflection of num 29
 num 31     : ["L-", "R+", "S+", "R+"], arguments --  c_c2scb(-x,y,-phi,am,b2)      -> time flip
 num 32     : ["R-", "L+", "S+", "L+"], arguments --  c_c2scb(-x,-y,phi,ap,b2)      -> time flip + reflection

'''

path_dict[25] = ["L+", "R-", "S-", "L-"]
path_dict[26] = ["R+", "L-", "S-", "R-"]
path_dict[27] = ["L-", "R+", "S+", "L+"]
path_dict[28] = ["R-", "L+", "S+", "R+"]

path_dict[29] = ["L+", "R-", "S-", "R-"]
path_dict[30] = ["R+", "L-", "S-", "L-"]
path_dict[31] = ["L-", "R+", "S+", "R+"]
path_dict[32] = ["R-", "L+", "S+", "L+"]

'''
 ########## ---------------------------------- C | C2 S C2 | C   Curves -------------------------------------------- 
 num 33     : ["L+", "R-", "S-", "L-", "R+"], arguments -- c_c2sc2_c(x,y,phi,ap,b2)       -> given 
 num 34     : ["R+", "L-", "S-", "R-", "L+"], arguments -- c_c2sc2_c(x,-y,-phi,am,b2)     -> reflection of num 33
 num 35     : ["L-", "R+", "S+", "L+", "R-"], arguments -- c_c2sc2_c(-x,y,-phi,am,b2)     -> time flip of num 33
 num 36     : ["R-", "L+", "S+", "R+", "L-"], arguments -- c_c2sc2_c(-x,-y,phi,ap,b2)      -> time flip of num 34     

'''

path_dict[33] = ["L+", "R-", "S-", "L-", "R+"]
path_dict[34] = ["R+", "L-", "S-", "R-", "L+"]
path_dict[35] = ["L-", "R+", "S+", "L+", "R-"]
path_dict[36] = ["R-", "L+", "S+", "R+", "L-"]

'''
 ########## ---------------------------------- C C | C    Curves --------------------------------------------------- 

  num 37    : ["L+", "R+", "L-"], arguments -- cc_c(x,y,phi,ap,b1)        -> given
  num 38    : ["R+", "L+", "R-"], arguments -- cc_c(x,-y,-phi,am,b1)      -> reflection of num 37
  num 39    : ["L-", "R-", "L+"], arguments -- cc_c(-x,y,-phi,am,b1)      -> time flip of num 37 
  num 40    : ["R-", "L-", "R+"], arguments -- cc_c(-x,-y,phi,ap,b1)      -> time flip + reflection

'''

path_dict[37] = ["L+", "R+", "L-"]
path_dict[38] = ["R+", "L+", "R-"]
path_dict[39] = ["L-", "R-", "L+"]
path_dict[40] = ["R-", "L-", "R+"]

'''
 ########## ----------------------------------  C S C2 | C Curves --------------------------------------------------
 num 41     : ["L+", "S+", "R+", "L-"], arguments -- csc2_ca(x,y,phi,ap,b1)         -> given
 num 42     : ["R+", "S+", "L+", "R-"], arguments -- csc2_ca(x,-y,-phi,am,b1)       -> reflection
 num 43     : ["L-", "S-", "R-", "L+"], arguments -- csc2_ca(-x,y,-phi,am,b1)       -> time flip
 num 44     : ["R-", "S-", "L-", "R+"], arguments -- csc2_ca(-x,-y,phi,ap,b1)       -> time flip + reflection

 num 45     : ["L+", "S+", "L+", "R-"], arguments -- csc2_cb(x,y,phi,ap,b2)         -> given
 num 46     : ["R+", "S+", "R+", "L-"], arguments -- csc2_cb(x,-y,-phi,am,b2)       -> reflection
 num 47     : ["L-", "S-", "L-", "R+"], arguments -- csc2_cb(-x,y,-phi,am,b2)       -> time flip
 num 48     : ["R-", "S-", "R-", "L+"], arguments -- csc2_cb(-x,-y,phi,ap,b2)       -> time flip and reflection


'''
path_dict[41] = ["L+", "S+", "R+", "L-"]
path_dict[42] = ["R+", "S+", "L+", "R-"]
path_dict[43] = ["L-", "S-", "R-", "L+"]
path_dict[44] = ["R-", "S-", "L-", "R+"]

path_dict[45] = ["L+", "S+", "L+", "R-"]
path_dict[46] = ["R+", "S+", "R+", "L-"]
path_dict[47] = ["L-", "S-", "L-", "R+"]
path_dict[48] = ["R-", "S-", "R-", "L+"]
