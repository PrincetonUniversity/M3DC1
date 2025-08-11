# Simmetrix to M3DC1 model/mesh writer (simToM3dc1).
The simToM3dc1 reads a Simmetrix model/mesh, writes PUMI model/mesh along
with the normal vector and curvature data on mesh vertices of selected loops.
Also writes a .txt file (Output file name\_modelInfo.txt) that contains the following information
  * Minimum of bounding box (x\_min, y\_min)
  * Maximum of bounding box (x\_max, y\_max)
  * Number of model edges on Innermost loop
  * model edges Ids on Innermost loop
  * Number of model edges on Outermost loop
  * model edges Ids on Outermost loop

* Takes in an input file that contains:
    * Simmetrix model (.smd)
    * Simmetrix mesh(.sms)  
    * Native Parasolid(.x\_t, .xmt\_txt, .x\_b, .xmt\_bin) (Not optional if model is created from Parasolid).
    * Output file name (optional)
    * Innermost model face number from Simmetrix model
    * Outermost model face number from Simmetrix model

* Outputs
  * PUMI model (.dmg)
  * PUMI mesh (.smb)
  * A .txt file (Output file name\_modelInfo.txt) that contains the following additional model information
    * Minimum of bounding box (x\_min, y\_min)
    * Maximum of bounding box (x\_max, y\_max)
    * Number of model edges on Innermost loop
    * model edges Ids on Innermost loop
    * Number of model edges on Outermost loop
    * model edges Ids on Outermost loop

* Usage:
    * /path to directory with simToM3dc1/simToM3dc1 inputFileName

