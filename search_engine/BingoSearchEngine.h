#pragma once

#include "indigo.h"
#include "bingo-nosql.h"

#include "IndigoMolecule.h"
#include "IndigoQueryMolecule.h"
#include "IndigoSession.h"
#include "IndigoSDFileIterator.h"
#include "IndigoException.h"
#include "BingoNoSQL.h"
#include "BingoResultIterator.h"

#include "SearchEngineInterface.h"

#include <cstring>

class BingoSearchEngine : public SearchEngineInterface {
public:
    BingoSearchEngine() = delete;
    explicit BingoSearchEngine(const indigo_cpp::IndigoSessionPtr &);

    ~BingoSearchEngine() override;

    /**
     * @param path If path ends with ".sdf", searchEngine will be created from .sdf file,
     * otherwise it will try to load from path dir
     */
    void build(const std::string &path) override;

    std::vector<indigo_cpp::IndigoMolecule> findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) override;

private:
    int _db = -1;
    indigo_cpp::IndigoSessionPtr _indigoSessionPtr;
};