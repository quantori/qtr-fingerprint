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
    inline explicit IdConverter(const std::filesystem::path &idToStringDirPath) {
        LOG(INFO) << "Start id converter loading";
        for (auto &filename: qtr::findFiles(idToStringDirPath, ".csv")) {
            qtr::IdToStringReader reader(filename);
            for (const auto &[dbId, outerId]: reader)
                _fromDbId[dbId] = {outerId, _libraryIds.size()};
            _libraryIds.emplace_back(filename.stem());
        }
        LOG(INFO) << "Finish id converter loading";
    }

    inline std::pair<const std::string &, const std::string &> fromDbId(uint64_t innerId) const {
        auto it = _fromDbId.find(innerId);
        assert(it != _fromDbId.end());
        const auto &[outerId, library] = it->second;
        return {outerId, _libraryIds[library]};
    }
};
