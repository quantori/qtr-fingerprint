#pragma once

#include <string>
#include <unordered_map>

#include "glog/logging.h"

namespace qtr {

    enum class DbType {
        OnDrive,
        InRamSmiles,
        InRamMolecules
    };

    static std::unordered_map<std::string, DbType> dbTypeMap = {
            {"in_ram_smiles",    DbType::InRamSmiles},
            {"in_ram_molecules", DbType::InRamMolecules},
            {"on_drive",         DbType::OnDrive}
    };

    inline DbType parseDbType(const std::string &str) {
        if (!dbTypeMap.contains(str)) {
            LOG(ERROR) << "Bad mode option value";
            exit(-1);
        }
        return dbTypeMap[str];
    }

    enum class RunModeType {
        Interactive,
        FromFile,
        Web
    };

    static std::unordered_map<std::string, RunModeType> runModeMap = {
            {"interactive", RunModeType::Interactive},
            {"from_file",   RunModeType::FromFile},
            {"web",         RunModeType::Web}
    };

    inline RunModeType parseRunMode(const std::string &str) {
        if (!runModeMap.contains(str)) {
            LOG(ERROR) << "Bad run mode option";
            exit(-1);
        }
        return runModeMap[str];
    }

}