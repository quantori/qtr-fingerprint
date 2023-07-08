
#include "ArgsBase.h"
#include "stdexcept"
#include "Utils.h"

using namespace std;

namespace qtr {
    ArgsBase::ArgsBase(int argc, char **argv) {
        absl::ParseCommandLine(argc, argv);
        parseDbType();
        checkDbType();
        parseDbName();
        checkDbName();
    }

    void ArgsBase::parseDbType() {
        string dbTypeStr = GET_FLAG(dbType);
        if (dbTypeStr == "qtr_ram")
            _dbType = DataBaseType::QtrRam;
        else if (dbTypeStr == "qtr_drive")
            _dbType = DataBaseType::QtrDrive;
        else if (dbTypeStr == "bingo_nosql")
            _dbType = DataBaseType::BingoNoSQL;
        else {
            _dbType = DataBaseType::BadType;
            LOG(ERROR) << "Bad db_type option value: \"" << dbTypeStr << "\"";
            exit(-1);
        }
        LOG(INFO) << OPTION_NAME(dbType) << ": " << dbTypeStr;
    }

    void ArgsBase::parseDbName() {
        _dbName = GET_FLAG(dbName);
        LOG(INFO) << OPTION_NAME(dbName) << ": " << _dbName;
    }


    void ArgsBase::checkDbName() {
        checkEmptyArgument(_dbName, "Please specify " OPTION_NAME(dbName) " option");
    }

    void ArgsBase::checkDbType() {
        if (_dbType == DataBaseType::BadType) {
            LOG(ERROR) << "Please specify " << OPTION_NAME(_dbType) << " option";
            exit(-1);
        }
    }


    ArgsBase::DataBaseType ArgsBase::dbType() const {
        if (_dbType == DataBaseType::BadType) {
            LOG(ERROR) << tryToGetUninitializedField << OPTION_NAME(dbType);
            exit(-1);
        }
        return _dbType;
    }

    std::string ArgsBase::dbName() const {
        if (_dbName.empty()) {
            LOG(ERROR) << tryToGetUninitializedField << OPTION_NAME(dbName);
            exit(-1);
        }
        return _dbName;
    }

    ArgsBase::DataBaseType ArgsBase::stringToDataBaseType(const string &s) {
        if (s == "qtr_ram")
            return DataBaseType::QtrRam;
        else if (s == "qtr_drive")
            return DataBaseType::QtrDrive;
        else if (s == "bingo_nosql")
            return DataBaseType::BingoNoSQL;
        else {
            return DataBaseType::BadType;
        }
    }


} // qtr