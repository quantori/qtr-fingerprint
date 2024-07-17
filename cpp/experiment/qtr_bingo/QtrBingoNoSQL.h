#pragma once

#include "BingoResultIterator.h"
#include "IndigoMolecule.h"
#include "IndigoQueryMolecule.h"
#include "IndigoSession.h"
#include "IndigoSimilarityMetric.h"
#include "BingoNoSQL.h"
#include "qtr-bingo-nosql.h"
#include "base_cpp/array.h"

class QtrBingoNoSQL {
public:
    using target_t = indigo_cpp::IndigoMolecule;
    using query_t = indigo_cpp::IndigoQueryMolecule;

    int insertRecord(const target_t& entity);

    static QtrBingoNoSQL createDatabaseFile(indigo_cpp::IndigoSessionPtr session, const std::string& path, const std::string& options = "");

    indigo_cpp::BingoResultIterator<target_t>
    searchSub(const query_t &query, const indigo::Array<byte> &fp, const std::string &options = "") const;

    indigo_cpp::BingoResultIterator<target_t> searchSub(const query_t &query, const std::string &options = "") const;

    indigo_cpp::IndigoSessionPtr session;

private:
    QtrBingoNoSQL(indigo_cpp::IndigoSessionPtr indigo, int e);

    int id;
};
