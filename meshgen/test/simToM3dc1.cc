#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <set>
#include <cmath>
#include <cassert>
#include "SimUtil.h"
#include "SimModel.h"
#include "SimAdvModel.h"
#include "SimParasolidKrnl.h"
#include <SimPartitionedMesh.h>
#include <PCU.h>
#include "gmi_sim.h"
#include "gmi_mesh.h"
#include <apfMDS.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>

#ifdef LICENSE
  #include <SimLicense.h>
  #ifdef STELLAR
    char simLic[128]="/home/PPPL/simmetrix/license/simmetrix.lic";
  #endif
  #ifdef PPPL
    char simLic[128]="/opt/hpc/software/Simmetrix/simmodsuite.lic";
  #endif
  #ifdef SDUMONT
    char simLic[128]="/scratch/ntm/software/Simmetrix/license/simmodsuite.lic";
  #endif
  #ifdef MIT
    char simLic[128]="/orcd/nese/psfc/001/software/simmetrix/RLMServer-14/server.lic";
  #endif
#else // scorec
  char simLic[128]="/net/common/meshSim/license/license.txt";
#endif

/******************************************************
 ********************* Structs ***********************/

// Struct to hold input data.
struct Inputs{
  std::string simModelName;  // input Simmetrix model
  std::string simMeshName;  // input Simmetrix mesh
  std::string outputName = "output";  // output file name
  std::string nativeModelName;  // native model
  int innerMostFace;  // innermost model face
  int outerMostFace;  // outermost model face
};
 
// A 2D vector struct holding dx, dy.
struct Vec{
  double dx = 0.0;
  double dy = 0.0;
};


// A 2D point struct holding x,y coordinates.
struct Point{
  double x;
  double y;
};

// Model loop struct containing all model vertices, and edges
// on it. 
struct Loop{
  std::vector <pGEdge> modelEdges; 
  std::vector <pGVertex> modelVertices;
};

/*****************************************************************************
 ************************** Function Definitions ****************************/

/*
 * Function to verify the inputs.
 * Returns the struct Inputs that contains all the inputs.
 */
Inputs verifyInputs(int argc, char** argv);

/*
 * Function to load Simmetrix model given model name and native parasolid model.
 * std::string& simModel (in): Input Simmetrix model name (.smd)
 * std::string& nativeModel (in): Input native parasolid model name. 
 * Returns Simmetrix pGModel.
 */
pGModel loadSimModel(const std::string& simModel,const std::string& nativeModel);

/*
 * Find the center of domain of pGModel. 
 * pGModel& simModel (in): Input Simmetrix model.
 * int coreFace (in)(optional): Plasma core model face number.
 * Returns center point (x,y).
 */
Point getDomainCenter(const pGModel& model, int coreFace = 0);

/*
 * Given the innermost model face, get the innermost model loop.
 * pGModel& simModel (in): Input Simmetrix model.
 * int& modelNumber (in): Innermost model face number.
 * Returns innermost loop.
 */
Loop getInnerMostLoop(const pGModel& model, const int& modelNumber);

/*
 * Given the outermost model face, get the outermost model loop.
 * pGModel& simModel (in): Input Simmetrix model.
 * int& modelNumber (in): Outermost model face number.
 * Returns outermost loop.
 */
Loop getOuterMostLoop(const pGModel& model, const int& modelNumber);

/*
 * Given a mesh and model loop, attach data (normal & curvature) on the mesh entities
 * classified on that model loop.
 * pMesh& mesh (in): Simmetrix mesh.
 * Loop loop (in): The geometric model loop on which data is desired.
 * Point center (in): Center of the domain.
 * pMeshDataId* dataId (in-out): Desired data id to be attached on mesh.
 */
void dataOnVerticesAtLoop(const pMesh& mesh, Loop loop, Point center, pMeshDataId* dataId);

/*
 * Function to get a vector of mesh vertices classified on geometric edge (ge).
 * pMesh m (in): Simmetrix mesh.
 * pGEdge ge (in): Geometric edge on the model. 
 * Returns a vector of mesh vertices classified on ge in mesh m.
 */
std::vector <pVertex> getMeshVerticesOnModelEdge(pMesh m, pGEdge ge);

/*
 * Given a model loop (edges already set), set an ordered list of model vertices
 * on the loop.
 * Loop& l (in-out): Loop on which ordered list of model vertices is desired.
 */
void setModelVerticesOnLoop(Loop& l);

/*
 * Given a mesh vertex v on model edge ge, this function returns normal vector on it.
 * pVertex v (in): Mesh vertex on which normal vector is desired.
 * pGEdge ge (in): Model edge on which mesh vertex v is classified on.
 * Loop gl (in): Model loop on which model edge is classified on.
 * Point center (in): Center of the domain. Needed to check direction of normal vectors.
 * Returns outward normal vector on the mesh vertex v. 
 */
Vec normalVecOnModelEdge(pVertex v, pGEdge ge, Loop gl, Point center);

/*
 * Given a mesh vertex v on model vertex gv, this function returns normal vector on it.
 * pVertex v (in): Mesh vertex on which normal vector is desired.
 * pGVertex gv (in): Model vertex on which mesh vertex v is classified on.
 * Loop gl (in): Model loop on which model vertex is classified on.
 * Point center (in): Center of the domain. Needed to check direction of normal vectors.
 * Returns outward normal vector on the mesh vertex v.
 */
Vec normalVecOnModelVertex(pVertex v, pGVertex gv, Loop gl, Point center);

/*
 * Given a mesh vertex v on model edge ge, this function returns curvature on it.
 * pVertex v (in): Mesh vertex on which curvature is desired.
 * pGEdge ge (in): Model edge on which mesh vertex v is classified on.
 * Returns curvature value.
 */
double curvatureOnMeshVertex(pVertex v, pGEdge ge);

/*
 * Given a closed loop l, this function finds the area enclosed by the loop.
 * Since we only use end points, it might not be accurate. but it works for
 * our case where we need to compare areas of concentric loops.
 * Returns area value.
 */
double areaOfPolygon(Loop l);

/*
 * Given a point pt, and Loop l, this function checks if the point lies on the area enclosed
 * by the loop. Use crossing method to test.
 * Point pt (in): Test point.
 * Loop l (in): Loop that defines the area enclosed under it.
 * Returns true if pt lies on the area enclosed by the loop.
 */
bool pointInsideLoop(Point pt, Loop l);

/*
 *  Given the name of Simmetrix model, and an output model name, write (.dmg) model to disk
 *  and gmi model in memory for further use.
 *  std::string& simModel (in): Input Simmetrix model name (.smd).
 *  std::string output (out): Output model name to be written on disk. 
 *  Returns gmi_model.
 */
gmi_model* writePumiModel(std::string simModel, std::string output);

/*
 * Given the name of Simmetrix Mesh, Simmetrix pGModel, PUMI model (gmi),and  an output mesh 
 * name, write (.smb) mesh to disk.
 * std::string simMesh (in): Input Simmetrix mesh  name (.sms).
 * pGModel simModel (in): Simmetrix model.
 * gmi_model* mdl (in): PUMI model.
 * std::string output (out): Output mesh name to be written on disk.
 */
void writePumiMesh(std::string simMesh, pMesh m, pMeshDataId dataId, 
                  pGModel simModel, gmi_model* mdl, std::string output);

/*
 * A function to write additional model information to a txt file
 * pGModel simModel (in): Simmetrix model.
 * Loop innerLoop (in): Inner model loop to be written to file.
 * Loop outerLoop (in): Outer model loop to be written to file.
 * std::string outModelInfo (out): txt file containing model info  to be written on disk.
 */
void writeModelInfo(pGModel model, Loop innerLoop, Loop outerLoop, std::string outModelInfo);

/*******************************************************************************************
 ************************************** Main Function *************************************/
int main(int argc, char** argv)
{
  // Step 0: Initialize the MPI, PCU environment
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();

  // Step 1: Verify the inputs
  Inputs in = verifyInputs(argc, argv);

  // Step 2: Find Simmetrix license 
  #ifdef LICENSE
    SimLicense_start("geomsim_core,geomsim_adv,meshsim_surface,meshsim_adv",simLic);
  #else
    // for SCOREC
    Sim_readLicenseFile(simLic);
  #endif

  // Step 3: Load Simmetrix model and mesh
  pGModel model = loadSimModel(in.simModelName, in.nativeModelName);
  std::cout << "Loading Simmetrix mesh (.sms) ...\n";
  pMesh mesh = M_load(in.simMeshName.c_str(), model, NULL);
  std::cout << "Loaded mesh successfully\n";

  // Step 4: Read data from Simmetrix model
  std::cout << "Setting innermost and outermost loops\n";
  Point center = getDomainCenter(model, in.innerMostFace);
  Loop innerLoop = getInnerMostLoop(model, in.innerMostFace);
  Loop outerLoop = getOuterMostLoop(model, in.outerMostFace);

  // Step 5: Find normal vectors & curvatures at mesh vertices 
  // of given loops
  std::cout << "Setting normal vector and curvature data to mesh vertices on defined loops\n";
  pMeshDataId norm_curv = MD_newMeshDataId("norm_curv_id");
  dataOnVerticesAtLoop(mesh, innerLoop, center, &norm_curv);
  dataOnVerticesAtLoop(mesh, outerLoop, center, &norm_curv);

  // Step 6: Write Pumi model (.dmg), mesh (.smb), modelInfo(.txt)
  std::string outModel = in.outputName + ".dmg";
  std::string outMesh = in.outputName + ".smb";
  std::string outModelInfo = in.outputName + "_modelInfo.txt";
  gmi_model* mdl = writePumiModel(in.simModelName, outModel); 
  writePumiMesh(in.simMeshName, mesh, norm_curv, model, mdl, outMesh);
  writeModelInfo(model, innerLoop, outerLoop, outModelInfo);

  // Step 7: Terminate the MPI, PCU environments.
  PCU_Comm_Free();
  MPI_Finalize();
}

// Function to check if a string has given extension
// or not.
bool hasExtension(std::string s, std::string ext) 
{
  if(s.substr(s.find_last_of(".") + 1) == ext) 
    return true;
  else 
    return false;
}

// Function to verify the inputs.
Inputs verifyInputs(int argc, char** argv)
{
  // Step 1: Verify we have the input file
  if (argc < 2)
  {
    std::cout << "Inputs Missing\n";
    std::cout << "Usage: ./simToM3dc1 input_filename\n";
    exit(1);
  }

  // Step 2: Read the input file
  std::string inputFileName = argv[1];
  std::ifstream inFile(inputFileName);
  if (!inFile.is_open())
  {
    std::cout << "Error opening the input file " <<  argv[1] << "\n";
    exit(1);
  }

  // Step 3: Read the inputs and populate the Inputs struct
  Inputs in;
  std::string token;

  bool simModelFound = false;
  bool simMeshFound = false;
  bool innerFaceFound = false;
  bool outerFaceFound = false;

  while (inFile >> token)
  {
    // Step 3.1: Verify we have correct Simmetrix Model in input file.
    if (token == "sim_model")
    {
      inFile >> in.simModelName;
      simModelFound = true;
      if (!hasExtension(in.simModelName, "smd"))
      {
        std::cout << "Simmetrix model with extension .smd could not be found\n";
        std::cout << "Check model name provided for sim_model in the input file\n";
        exit(1);
      }
    }

    // Step 3.2: Verify we have correct Simmetrix Mesh in input file.
    if (token == "sim_mesh")
    {
      inFile >> in.simMeshName;
      simMeshFound = true;
      if (!hasExtension(in.simMeshName, "sms"))
      {
        std::cout << "Simmetrix mesh with extension .sms could not be found\n";
        std::cout << "Check mesh name provided for sim_mesh in the input file\n";
        exit(1);
      }
    }
    
    // Step 3.3: Read output file name
    if (token == "out_file")
      inFile >> in.outputName;

    // Step 3.4: Read native model file (parasold)
    if (token == "native_model")
      inFile >> in.nativeModelName;

    // Step 3.5: Read inner most model face number (simmetrix model face)
    if (token == "inner_face")
    {
      inFile >> in.innerMostFace;
      innerFaceFound = true;
    }
    
    // Step 3.6: Read outer most model face number (simmetrix model face)
    if (token == "outer_face")
    {
      inFile >> in.outerMostFace;
      outerFaceFound = true;
    }
  }
  inFile.close();

  // Step 4: Check if some required parameter is missing or not.
  if (!simModelFound)  // Check Simmetrix Model
  {
    std::cout << "Parameter sim_model missing from input file\n";
    exit(1);
  }
 
  if (!simMeshFound)  // Check Simmetrix Mesh
  {
    std::cout << "Parameter sim_mesh missing from input file\n";
    exit(1);
  }

  if (!innerFaceFound)  // Check Innermost Face Number
  {
    std::cout << "Parameter inner_face missing from input file\n";
    exit(1);
  }

  if (!outerFaceFound)  // Check Outermost Face Number
  {
    std::cout << "Parameter inner_face missing from input file\n";
    exit(1);
  }

  // Step 5: Return Inputs (in) after validating inputs
  return in;
}

// Function to load Simmetrix model given model name and native parasolid model.
pGModel loadSimModel(const std::string& simModel, const std::string& nativeModel)
{
  std::cout << "Loading Simmetrix Model (.smd) ...\n";
  // Step 1: Read Simmetrix model name, and see if native parasolid model
  // is provided or not.
  const char* simMdl = simModel.c_str();
  const char* nativeMdl;
  if (!nativeModel.empty())
    nativeMdl = nativeModel.c_str();
  else
    nativeMdl = nullptr;    

  // Step 2: Initialize Simmetrix objects.
  pProgress prog = Progress_new();
  Progress_setDefaultCallback(prog);
  SimParasolid_start(1);

  // Step 3: If native model is provided, read it and verify that extension is 
  // correct for parasolid models. 
  pNativeModel nModel = NULL;
  if (nativeMdl != nullptr)
  {
    int fileFormat = 0;  // Parasolid native file format
    
    // Step 3.1: Decide the file format based on the extension of given model.
    if (hasExtension(nativeModel,"x_t") || hasExtension(nativeModel,"xmt_txt"))
      fileFormat = 0;
    if (hasExtension(nativeModel,"x_b") || hasExtension(nativeModel,"xmt_bin"))
      fileFormat = 1;

    // Step 3.2: Using the file format (0 or 1), and native model name, load
    // native model.
    nModel = ParasolidNM_createFromFile(nativeMdl,fileFormat);
  }
  else  // If not native model, issue a warning.
  {
    std::cout << "Warning: If Simmetrix model was created from Parasolid model, make sure\n"
                  "to provide native parasolid model (.x_t, .xmt_txt, .x_b, .xmt_bin) \n";
  }  

  // Step 4: Load the Simmetrix geometric model. Should work with and without native model 
  // provided. If model was created without native model, then native model is not required.
  // If created with native model, code will throw errors in query functions, if native 
  // model isn't provided.
  pGModel geomModel = GM_load(simMdl, nModel, prog);
  std::cout << "Loaded model Successfully\n";
  
  // Step 5: Return the loaded Simmetrix model.
  return geomModel;  
}

// Function to find approximate center of domain. If core face is given,
// its centroid is used. otherwise, center of bounding box.
Point getDomainCenter(const pGModel& model, int coreFace)
{
  // Step 1: Find bounding box of the model.
  double min[3], max[3];
  GM_bounds(model, &min[0], &max[0]);

  // Step 2: Find center of bounding box.
  Point center;
  center.x = (min[0]+max[0])/2;
  center.y = (min[2]+max[2])/2;

  // Step 3: If innermost face is given, use its centroid instead.
  if (coreFace > 0)
  {
    double centroid[3];

    // Step 3.1: Get the model face from tag.
    pGFace gf = GM_faceByTag(model, coreFace);

    // Step 3.2: Find centroid and set it to return point.
    GF_centroid(gf, 1.0, centroid);
    center.x = centroid[0];
    center.y = centroid[1];
  }
  return center;
}

// Given the innermost model face, get the innermost model loop.
Loop getInnerMostLoop(const pGModel& model, const int& modelNumber)
{
  Loop inner;  // return loop.

  // Step 1: Get the model face using its tag, and find number of model loops on it.
  pGFace gf = GM_faceByTag(model, modelNumber);
  int numLoops = GF_numLoops(gf);
  assert (numLoops == 1);

  // Step 2: Get the loop iterator from model face (gf) and set the model edges to
  // the return loop.
  int side = 1;
  pGFaceUse fu = GF_use(gf,side);
  pGLUIter loopIter = GFU_loopIter(fu);
  while (pGLoopUse loopUse = GLUIter_next(loopIter))
  {
    // Step 2.1: Iterate over model edges on each loop and set these edges
    // to the return loop.
    pGEUIter edgesOnLoop = GLU_edgeUseIter(loopUse);
    while (pGEdgeUse edgeUse = GEUIter_next(edgesOnLoop))
    {
      // Step 2.1.1: Read each model edge and push it to the model edges 
      // vector in the definition of loop.
      pGEdge ge =  GEU_edge(edgeUse);
      inner.modelEdges.push_back(ge);
    }
  }
  GLUIter_delete(loopIter);

  // Step 3: Based on model edges defined, set the ordered list of model vertices
  // to each loop.
  setModelVerticesOnLoop(inner);

  return inner;
}

// Given the outermost model face, get the outermost model loop. 
Loop getOuterMostLoop(const pGModel& model, const int& modelNumber)
{
  Loop outer;  // return loop
  
  // Step 1: Get the model face using its tag.
  pGFace gf = GM_faceByTag(model, modelNumber);

  // Step 2: Get the loop iterator from model face (gf), and set the model edges and vertices
  // to it. Decide which loop is the outermost of all the loops on the face.
  int side = 1;
  pGFaceUse fu = GF_use(gf,side);

  // Step 2.1: Iterate over the loops on the model face.
  pGLUIter loopIter = GFU_loopIter(fu); 
  double area = 0.0;  
  while (pGLoopUse loopUse = GLUIter_next(loopIter))
  {
    // Step 2.2: Use temp loop until find the outermost one. Get the model edges
    // on the loop.
    Loop temp;
    pGEUIter edgesOnLoop = GLU_edgeUseIter(loopUse);
    while (pGEdgeUse edgeUse = GEUIter_next(edgesOnLoop))
    {
      // Step 2.2.1: Read each model edge and push it to the model edges 
      // vector in the definition of loop.
      pGEdge ge =  GEU_edge(edgeUse);
      temp.modelEdges.push_back(ge);
    }
    
    // Step 2.3: Based on model edges defined, set the ordered list of model vertices to the 
    // temp loop.
    setModelVerticesOnLoop(temp);

    // Step 2.4: Find the loop with maximum under under curve and set it to outer loop.
    double areaPolygon = areaOfPolygon(temp);
    if (fabs(areaPolygon) > area)
    { 
      outer = temp;
      area = fabs(areaPolygon);
    }
  }
  GLUIter_delete(loopIter);


  return outer;
}

// Function to return a vector of mesh vertices classified on geometric edge (ge).
std::vector <pVertex> getMeshVerticesOnModelEdge(pMesh m, pGEdge ge)
{
  std::vector <pVertex> vertices;

  // Simmetrix mesh vertex iterator on model edge ge.
  VIter vertexIter = M_classifiedVertexIter(m, ge, 1);
  while (pVertex v  = VIter_next(vertexIter))
    vertices.push_back(v);
  VIter_delete(vertexIter);
  
  return vertices;
}

// Given a model loop (edges already set), set an ordered list of model vertices
// on the loop.
void setModelVerticesOnLoop(Loop& l)
{
  std::vector <pGEdge> modelEdges = l.modelEdges;
  
  // Step 1: Find the model vertices on model edge 1 and 2 to figure out 
  // correct order of the vertices.
  pGVertex gv0 = GE_vertex(modelEdges[0],0);
  pGVertex gv1 = GE_vertex(modelEdges[0],1);
  pGVertex gv2 = GE_vertex(modelEdges[1],0); 
  pGVertex gv3 = GE_vertex(modelEdges[1],1);

  // Step 2: If gv0 is the common vertex between two edges, the vector will 
  // start with gv1 and then gv0  and so on: Otherwise, gv0 is the first 
  // vertex and gv1 is next and so on. 
  if (gv0 == gv2 || gv0 == gv3)
  {
    l.modelVertices.push_back(gv1);
    l.modelVertices.push_back(gv0);
  }
  else
  {
    l.modelVertices.push_back(gv0);
    l.modelVertices.push_back(gv1);
  }

  // Step 3: Once we determine the order of first two vertices, from there
  // its only looking for the common vertex between two model edges. This is
  // done by comparing the two vertices on the model edge with the last saved
  // vertex in the vertices vector.
  for (int i = 1; i < modelEdges.size(); i++)
  {
    gv0 = GE_vertex(modelEdges[i],0);
    gv1 = GE_vertex(modelEdges[i],1);
    
    pGVertex gvLast = l.modelVertices[l.modelVertices.size()-1];
    if(gv0 == gvLast)
      l.modelVertices.push_back(gv1);
    else
      l.modelVertices.push_back(gv0);
  }
}

// Given a mesh and model loop, attach data (normal & curvature) on the mesh entities classified
// on the model loop.
void dataOnVerticesAtLoop(const pMesh& mesh, Loop loop, Point center, pMeshDataId* dataId)
{
  // Step 1: Iterate over the model edges in the loop.
  std::vector <pGEdge> modelEdges = loop.modelEdges;
  for (int i = 0; i < modelEdges.size(); i++)
  {
    // Step 2: For each model edge, get the mesh vertices classified on it.
    pGEdge ge = modelEdges[i];
    std::vector <pVertex> meshVertices = getMeshVerticesOnModelEdge(mesh, ge);

    // Step 3: Iterate over mesh vertices.
    for (int j = 0; j < meshVertices.size(); j++)
    {
      // Step 4: Depending on what model entity mesh vertex is classified on, find the 
      // normal vector.
      Vec norm;
      pVertex mV = meshVertices[j];
      gType entType = V_whatInType(mV);
      if (entType == 0)  // mesh vertex classified on model vertex
      {
        pGEntity gv = V_whatIn(mV);
        norm = normalVecOnModelVertex(mV, pGVertex(gv), loop, center);
      } 
      else if (entType == 1)  // mesh vertex classified on model edge
        norm = normalVecOnModelEdge(mV, ge, loop, center);

      // Step 5: Find the curvature on the mesh vertices.
      double curvature = curvatureOnMeshVertex(mV, ge);

      // Step 6: Attach data (normal vector & curvature)  to the mesh vertex.
      if (!EN_getDataPtr(mV, *dataId, NULL))
      {
        double *vData = new double[3];
        vData[0] = norm.dx;
        vData[1] = norm.dy;
        vData[2] = curvature;
        EN_attachDataPtr((pEntity)mV, *dataId, (void*)vData);
      }  // end if loop
    }  // end inner for loop (mesh vertices)  
  }  // end outer for loop (model edges)
}

// Given a mesh vertex v on model edge ge, this function returns normal vector on it.
Vec normalVecOnModelEdge(pVertex v, pGEdge ge, Loop gl, Point center)
{
  Vec norm;  // return vector (normal)

  // Step 1: Find coordinates of vertex v, and parametric range of model edge ge.
  pPoint pt = V_point(v);
  double parR[2];
  GE_parRange(ge, &parR[0], &parR[1]);
  
  // Step 2: Find the parametric value of mesh vertex depending on the correct 
  // classification type and assert that its within correct range. 
  double t = 0.0;  // parametric value
  gType entType = V_whatInType(v);
  if (entType == 0)  // classified on model vertex
  {
    pGEntity gv = V_whatIn(v);    
    t = GE_vertexReparam(ge, pGVertex(gv));
  }  
  else if (entType == 1)  // classified on modele edge
    t = P_param1(pt);

  assert(t >= parR[0] && t <= parR[1]);

  // Step 3: Find the tangent vector on parametric location t.
  double dt[3];
  GE_firstDerivative(ge, t, &dt[0]);

  // Step 4: x & y component of above derivate, and its magnitude.
  double dx = dt[0];
  double dy = dt[2];
  double len = sqrt(dx*dx+dy*dy);

  // Step 5: Find the vector from the vertex v to center point.
  double vec[2] = {center.x - P_x(pt), center.y - P_z(pt)}; 
  
  // Step 6: Find the normal vector (temporaray since direction needs to be adjusted)
  double normTemp[2];
  normTemp[0] = dy/len;
  normTemp[1] = -1.0*dx/len;
  
  // Step 7: Find the direction of normal vector (inward or outward)
  double dir = normTemp[0]*vec[0] + normTemp[1]*vec[1];

  // Step 8: Based on direction, find the final normal vector
  if (dir > 0) 
  {
    norm.dx = -normTemp[0];
    norm.dy = -normTemp[1];
  }
  else
  {
    norm.dx = normTemp[0];
    norm.dy = normTemp[1];	
  }  

  // Step 8: In complicated shapes, direction from centroid, and orientation-tangent might not work
  // perfectly to find if normal in inward or outward. Add this extra check at the end.
  double offset = 0.5;
  Point ptTest;
  ptTest.x = P_x(pt) + offset*norm.dx;
  ptTest.y = P_z(pt) + offset*norm.dy;
  
 if (pointInsideLoop(ptTest, gl))
  {
    norm.dx = -norm.dx;
    norm.dy = -norm.dy;
  }

  return norm;
}

// Given a mesh vertex v on model vertex gv, this function returns normal vector on it. 
Vec normalVecOnModelVertex(pVertex v, pGVertex gv, Loop gl, Point center)
{
  // Step 1: Find the model edges adjacent to model vertex gv, and also lie on the
  // loop gl.
  std::vector <pGEdge> edgesOnLoop;
  pPList edges = GV_edges(gv);
  for (int i = 0; i < PList_size(edges); i++)
  {
    pGEdge ge = static_cast<pGEdge>(PList_item(edges,i));
    for (int j = 0; j < gl.modelEdges.size(); j++)
    {
      if (ge == gl.modelEdges[j])
        edgesOnLoop.push_back(ge);
    }
  }
  PList_delete(edges);
  assert(edgesOnLoop.size() == 2);

  // Step 2: Find normal vector of the mesh vertex with respect to both the model edges
  // intersecting at gv.
  Vec norm0 = normalVecOnModelEdge(v, edgesOnLoop[0], gl, center);
  Vec norm1 = normalVecOnModelEdge(v, edgesOnLoop[1], gl, center);

  // Step 3: Average both the normal vectors for a final value.
  double dx = norm0.dx + norm1.dx;
  double dy = norm0.dy + norm1.dy;

  double len = sqrt(dx*dx+dy*dy);

  // Step 4: Unit vector of averaged normal vector.
  Vec norm;
  norm.dx = dx/len;
  norm.dy = dy/len;

  return norm;
}


double curvatureOnMeshVertex(pVertex v, pGEdge ge)
{
  // Step 1: Find the parametric value depending on the correct 
  // classification type and make sure its in the correct range. 
  double parR[2];
  GE_parRange(ge, &parR[0], &parR[1]);

  double t = 0.0;  // parametric value
  gType entType = V_whatInType(v);
  if (entType == 0)  // classified on model vertex
  {
    pGEntity gv = V_whatIn(v);    
    t = GE_vertexReparam(ge, pGVertex(gv));
  }
  else if (entType == 1)  // classified on model edge
  {
    pPoint pt = V_point(v);
    t = P_param1(pt);
  }
  assert(t >= parR[0] && t <= parR[1]);
 
  // Step 2: Find curvature at the given parametric value (t) at edge ge.
  double curvature;
  GE_curvature(ge, t, &curvature, NULL);
 
  return curvature;
}

// Return signed area of a polygon
double areaOfPolygon(Loop l)
{
  
  double area = 0.0;
  double pt0[3], pt1[3];
  pGVertex v0,v1;

  // Iterate over each segment between two model vertices, and find
  // area under curve.
  for (int i = 0; i < l.modelVertices.size()-1; i++)
  {
    v0 = l.modelVertices[i];
    v1 = l.modelVertices[i+1]; 
    
    GV_point(v0, pt0);
    GV_point(v1, pt1);

    area += (pt1[0] - pt0[0])*(pt1[2]+pt0[2]);
  }

  return area;
}

// Given a point pt, and Loop l, this function checks if the point lies on the area enclosed
// by the loop. Use ray crossing method to test.
bool pointInsideLoop(Point pt, Loop l)
{
  int crossings = 0;

  double pt0[3], pt1[3];
  pGVertex v0,v1;

  // Step 1: Iterate over model vertices (we need to look at a segment between two 
  // vertices in each loop)
  for (int i = 0; i < l.modelVertices.size()-1; i++)
  {
    pGVertex gv0 = l.modelVertices[i];
    pGVertex gv1 = l.modelVertices[i+1];

    GV_point(gv0, pt0);
    GV_point(gv1, pt1);

    // Step 2: If both points are on the opposide side of horizontal line drawn the 
    // test point pt, then find the intersection between that horizontal line and
    // line segment between (Pt[i]-Pt[i+1])
    if ((pt0[2] > pt.y) != (pt1[2] > pt.y))
    {
      double intersection = (pt1[0] - pt0[0])*(pt.y - pt0[2]) / (pt1[2] - pt0[2]) + pt0[0];
      if (pt.x < intersection)  // Check if test point is on right of intersection
        crossings++;
    }
  }
  
  // Step 3: If crossing i odd, point is inside.
  if (crossings % 2 == 1)
    return true;

  return false; 
}

// Given the name of Simmetrix model, and an output model name, write (.dmg) model 
// to disk and gmi model in memory for further use.
gmi_model* writePumiModel(std::string simModel, std::string output)
{
  std::cout << "Writing PUMI model (.dmg) ... \n";

  // Step 1: Initialize the needed objects.
  gmi_sim_start();
  gmi_register_mesh();
  gmi_register_sim();

  // Step 2: Get names written in strings as chars
  char simMdl[256], out[256]; 
  std::strncpy(simMdl, simModel.c_str(), sizeof(simMdl) - 1);
  std::strncpy(out, output.c_str(), sizeof(out) - 1);
 
  // Step 3: Write the model files
  gmi_model* mdl = gmi_load(simMdl);
  gmi_write_dmg(mdl, out);
  
  std::cout << "PUMI model written: " << output << "\n";
  
  gmi_sim_stop();
  return mdl;
}

// Write mesh data provided in pMeshDataId to the mesh vertices.
void writeDataOnPumiMesh(apf::Mesh2* mesh, pMesh simMesh, pMeshDataId dataId)
{
  // Step 1: Create a tag to write data. In our case, data is of type double 
  // and has three components (two normal vector and one curvature)
  pMeshDataId norm_curv = dataId;
  apf::MeshTag* norm_curv_tag = mesh->createDoubleTag("norm_curv", 3);

  // Step 2: Iterate over Simmetrix mesh and copy data written on Simmetrix mesh
  // vertices to corresponding PUMI mesh vertices.
  VIter vIter = M_vertexIter(simMesh);
  while (pVertex v =  VIter_next(vIter))
  {
    int id = EN_id(v);
    // Step 2.1: If data don't exist on vertex, do not take any action.
    // Else, set it to PUMI mesh vertex.
    if (EN_getDataPtr((pEntity)v, norm_curv, NULL))
    {
      // Step 2.1.1: Fetch the data from Simmetrix vertex.
      double* normCurvData;
      EN_getDataPtr((pEntity)v, norm_curv, (void**)&normCurvData);

      // Step 2.1.2: Attach the data to corresponding PUMI vertex.
      apf::MeshEntity* mV = getMdsEntity(mesh, 0, id);
      mesh->setDoubleTag(mV, norm_curv_tag, &normCurvData[0]);
    }
  }
  VIter_delete(vIter);  
}

// Given the name of Simmetrix Mesh, Simmetrix pGModel, PUMI model (gmi),
// and  an output mesh name, write (.smb) mesh to disk.
void writePumiMesh(std::string simMesh, pMesh m, pMeshDataId dataId, 
                  pGModel simModel, gmi_model* mdl, std::string output)
{
  std::cout << "Writing PUMI mesh (.smb) ... \n";

  // Step 1: Initialize the needed objects
  SimPartitionedMesh_start(NULL,NULL);
  gmi_sim_start();
  gmi_register_sim();

  // Step 2: Get names written in strings as chars
  char simMsh[256], out[256];
  std::strncpy(simMsh, simMesh.c_str(), sizeof(simMsh) - 1);
  std::strncpy(out, output.c_str(), sizeof(out) - 1);

  // Step 3: Write Simmetrix partitioned mesh and PUMIapf mesh from it.
  pParMesh simPMesh = PM_load(simMsh, simModel, NULL);
  apf::Mesh* simApfMesh = apf::createMesh(simPMesh);

  // Step 4: Write mesh data on PUMI mesh.
  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh); 
  writeDataOnPumiMesh(mesh, m, dataId);
  
  // Step 5: Write mesh on disk.
  mesh->writeNative(out);
  
  std::cout << "PUMI mesh written: " << output << "\n";
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  SimPartitionedMesh_stop();
}

// A function to write additional model information to a txt file.
// Write bounding box information and inner/outer loop information.
// Can add more information as we move forward on this.
void writeModelInfo(pGModel model, Loop innerLoop, Loop outerLoop, std::string outModelInfo)
{
  std::cout << "Writing additional model info (.txt) ... \n";  

  // Step 1: The file to write the additional model information.
  FILE* fp = fopen(outModelInfo.c_str(), "w");

  // Step 2: Find the bounding box the model and write it to the file:
  double min[3], max[3];
  GM_bounds(model, min, max);
  fprintf(fp, "%lf %lf\n", min[0], min[2]);
  fprintf(fp, "%lf %lf\n", max[0], max[2]);

  // Step 3: Write Inner Loop Info to the file.
  int numEdges = innerLoop.modelEdges.size();
  fprintf(fp, "%d\n" , numEdges);
  for (int i = 0 ; i < numEdges; i++)
    fprintf(fp, "%d ", GEN_tag(innerLoop.modelEdges[i]));

  fprintf(fp, "\n");
  // Step 3: Write Inner Loop Info to the file.
  numEdges = outerLoop.modelEdges.size();
  fprintf(fp, "%d\n" , numEdges);
  for (int i = 0 ; i < numEdges; i++)
    fprintf(fp, "%d ", GEN_tag(outerLoop.modelEdges[i]));
  
  std::cout << "Model info written: " << outModelInfo << "\n";
  fclose(fp);
}
