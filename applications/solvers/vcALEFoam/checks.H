//check if the mesh is tetrahedral
// get faces of cell 0
const labelList& face = mesh.cells()[0];

if ( face.size() == 4 ){
    //Tetrahedron
} else {
    FatalErrorIn("checks.H") << "The mesh is not tetrahedral." << abort(FatalError);
}
// return 0



//const word& solidModel__ = mechanicalProperties.subDict("solidModel").lookup("solidModel");
//
//
//
//
// check if FSI, if so: creates the dummy u
bool isFSI {false};
forAll(mesh.boundary(), patch) {
  const fvPatch& patx = mesh.boundary()[patch];
  if (patx.name().find("interface") != string::npos) { isFSI = true; } 
}


  pointVectorField uD(
        IOobject("uD", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), pMesh
  );

if (isFSI) {
  Info << "This is a FSI problem." << nl;
  uD = u;
}
