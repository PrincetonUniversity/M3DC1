# Simmetrix to M3DC1 model/mesh writer (simToM3dc1).
The simToM3dc1 reads a Simmetrix model/mesh, writes PUMI model/mesh alon
with the normal vector and curvature data on mesh vertices of selected loops.

* Takes in an input file that contains:
    * Simmetrix model (.smd)
    * Simmetrix mesh(.sms)  
    * Native Parasolid(.x\_t, .xmt\_txt, .x\_b, .xmt\_bin) (Not optional if model is created from Parasolid).
    * Output file name (optional)
    * Innermost model face number from Simmetrix model
    * Outermost model face number from Simmetrix model

* Usage:
    * /path to directory with simToM3dc1/simToM3dc1 inputFileName
