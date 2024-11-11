#pragma once

#include "GlobalIndigoSession.h"
#include "IndigoSubstructureMatcher.h"

inline bool isSubstructure(const indigo_cpp::IndigoQueryMolecule &queryMolecule,
                    const indigo_cpp::IndigoMolecule &candidateMol) {
    auto matcher = globalIndigoSession->substructureMatcher(candidateMol);
    return bool(indigoMatch(matcher.id(), queryMolecule.id()));
}