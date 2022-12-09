#pragma once

#include <unordered_map>
#include <string>
#include <vector>
#include <filesystem>

#include "Utils.h"
#include "id_to_string_io/IdToStringReader.h"

class IdConverter {
    std::vector<std::string> _libraryIds;
    std::unordered_map<uint64_t, std::pair<std::string, size_t>> _fromDbId; // from our db id to a pair of id and index in _libraryIds

public:
    inline explicit IdConverter(std::filesystem::path &idToStringDirPath) {
        for (auto &filename: qtr::findFiles(idToStringDirPath, ".csv")) {
            qtr::IdToStringReader reader(filename);
            for (const auto &[dbId, outerId]: reader)
                _fromDbId[dbId] = {outerId, _libraryIds.size()};
            _libraryIds.emplace_back(filename.stem());
        }
    }

    inline std::pair<std::string &, std::string &> fromDbId(uint64_t dbId) {
        return {_fromDbId.at(dbId).first, _libraryIds[_fromDbId.at(dbId).second]};
    }
};
