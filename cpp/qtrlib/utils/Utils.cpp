#include "Utils.h"

namespace qtr {

    bool endsWith(const std::string &a, const std::string &b) {
        return a.size() >= b.size() &&
               (0 == a.compare(a.size() - b.size(), b.size(), b));
    }

    int chexToInt(char letter) {
        return letter >= '0' && letter <= '9' ? letter - '0' : letter - 'a' + 10;
    }

    std::vector<std::filesystem::path> findFiles(const std::filesystem::path &pathToDir, const std::string &extension) {
        std::vector<std::filesystem::path> sdfFiles;
        for (const auto &entry: std::filesystem::recursive_directory_iterator(pathToDir)) {
            if (entry.path().extension() == extension) {
    //            LOG(INFO) << entry.path().string() << std::endl;
                sdfFiles.push_back(entry.path());
            }
        }
        return sdfFiles;
    }

    template <>
    void emptyArgument<uint64_t>(const uint64_t& argument, const std::string &message) {
        if (argument == 0) {
            LOG(ERROR) << message;
            exit(-1);
        }
    }

} // namespace qtr