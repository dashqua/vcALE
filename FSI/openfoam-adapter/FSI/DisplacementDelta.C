#include "DisplacementDelta.H"

using namespace Foam;

preciceAdapter::FSI::DisplacementDelta::DisplacementDelta
(
    const Foam::fvMesh& mesh,
    const std::string nameDPointDisplacement
)
:
DpointDisplacement_(
    const_cast<pointVectorField*>
    (
        &mesh.lookupObject<pointVectorField>(nameDPointDisplacement)
    )
)
{
    dataType_ = vector;
}

void preciceAdapter::FSI::DisplacementDelta::write(double * buffer, bool meshConnectivity, const unsigned int dim)
{
    // For every element in the buffer
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the displacement on the patch
        fixedValuePointPatchVectorField& DpointDisplacementFluidPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                DpointDisplacement_->boundaryFieldRef()[patchID]
            );

        // For every cell of the patch,
	// write displacements to precice buffer
        forAll(DpointDisplacement_->boundaryFieldRef()[patchID], i)
        {
            buffer[bufferIndex++] = DpointDisplacementFluidPatch[i][0];
            buffer[bufferIndex++] = DpointDisplacementFluidPatch[i][1];
            if(dim==3)
                buffer[bufferIndex++] = DpointDisplacementFluidPatch[i][2];
        }
    }
}

// return the displacement to use later in the velocity?
void preciceAdapter::FSI::DisplacementDelta::read(double * buffer, const unsigned int dim)
{
    // For every element in the buffer
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the displacement on the patch
        fixedValuePointPatchVectorField& pointDisplacementFluidPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                DpointDisplacement_->boundaryFieldRef()[patchID]
            );

        // For every cell of the patch
        forAll(DpointDisplacement_->boundaryFieldRef()[patchID], i)
        {
            // Add the received delta to the actual displacement
            pointDisplacementFluidPatch[i][0] = buffer[bufferIndex++];
            pointDisplacementFluidPatch[i][1] = buffer[bufferIndex++];
            if(dim==3)
                pointDisplacementFluidPatch[i][2] = buffer[bufferIndex++];
        }
    }
}
