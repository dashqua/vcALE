//check if the mesh is tetrahedral
// get faces of cell 0
const labelList& face = mesh.cells()[0];

if ( face.size() == 4 ){
    //Tetrahedron
} else {
    FatalErrorIn("checks.H") << "The mesh is not tetrahedral." << abort(FatalError);
}
// return 0



const word& solidModel__ = mechanicalProperties.subDict("solidModel").lookup("solidModel");
