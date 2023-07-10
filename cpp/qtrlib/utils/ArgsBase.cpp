
#include "ArgsBase.h"
#include "stdexcept"
#include "Utils.h"

using namespace std;

namespace qtr {
    ArgsBase::ArgsBase(int argc, char **argv) {
        absl::ParseCommandLine(argc, argv);
    }


    ArgsBase::DataBaseType ArgsBase::strToDataBaseType(const std::string &s) {
        auto it = _strToDataBaseType.find(s);
        if (it == _strToDataBaseType.end())
            return DataBaseType::BadType;
        return it->second;
    }

    std::vector<std::filesystem::path> ArgsBase::vecStrToVecPath(const vector <std::string> &v) {
        std::vector<std::filesystem::path> result;
        std::copy(v.begin(), v.end(), std::back_inserter(result));
        return result;
    }

} // qtr