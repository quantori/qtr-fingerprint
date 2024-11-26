#pragma once

#include "GlobalIndigoSession.h"
#include "IndigoSubstructureMatcher.h"
#include "Profiling.h"

inline bool isSubstructure(indigo::QueryMolecule &queryMolecule,
                           indigo::Molecule &candidateMol) {
    ProfileScope("isSubstructure");
    // TODO: maybe delete aromatize?
    candidateMol.aromatize(AromaticityOptions());
    MoleculeSubstructureMatcher msm(candidateMol);
    msm.setQuery(queryMolecule);
    return msm.find();
}