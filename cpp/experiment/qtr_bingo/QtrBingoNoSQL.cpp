#include "QtrBingoNoSQL.h"

#include <bingo-nosql.h>
#include "IndigoException.h"
#include "IndigoIterator.h"
#include "qtr-bingo-nosql.h"

namespace {
    template<typename target_t, typename query_t>
    constexpr const char *getDBTypeString() {
        if (std::is_same<target_t, indigo_cpp::IndigoMolecule>::value) {
            static_assert(std::is_same<query_t, indigo_cpp::IndigoQueryMolecule>::value, "");
            return "molecule";
        }
        //    else if (std::is_same<target_t, IndigoReaction>::value)
        //    {
        //        static_assert(std::is_same<query_t, IndigoQueryReaction>::value, "");
        //        return "reaction";
        //    }
        throw indigo_cpp::IndigoException("Unknown DB type");
    }
}

indigo_cpp::BingoResultIterator<QtrBingoNoSQL::target_t>
QtrBingoNoSQL::searchSub(const QtrBingoNoSQL::query_t &query, const std::string &options) const {
    session->setSessionId();
    return {session->_checkResult(bingoSearchSub(id, query.id(), options.c_str())), session};
}

indigo_cpp::BingoResultIterator<QtrBingoNoSQL::target_t>
QtrBingoNoSQL::searchSub(const QtrBingoNoSQL::query_t &query, const indigo::Array<byte> &fp,
                         const std::string &options) const {
    session->setSessionId();
    return {session->_checkResult(qtrBingoSearchSub(id, query.id(), options.c_str(), fp)), session};
}

int QtrBingoNoSQL::insertRecord(const QtrBingoNoSQL::target_t &entity) {
    session->setSessionId();
    return session->_checkResult(bingoInsertRecordObj(id, entity.id()));
}

QtrBingoNoSQL QtrBingoNoSQL::createDatabaseFile(indigo_cpp::IndigoSessionPtr session, const std::string &path,
                                                const std::string &options) {
    session->setSessionId();
    int id = session->_checkResult(
            bingoCreateDatabaseFile(path.c_str(), getDBTypeString<target_t, query_t>(), options.c_str()));
    return {std::move(session), id};
}

QtrBingoNoSQL::QtrBingoNoSQL(indigo_cpp::IndigoSessionPtr indigo, int e) : session(std::move(indigo)), id(e) {

}

