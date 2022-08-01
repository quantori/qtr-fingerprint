#pragma once

#include <functional>
#include <string>
#include <vector>

namespace qtr {

    using choose_func_t = std::function<std::vector<int>(const std::string &)>;

    class ColumnsChooser {
    public:
        ColumnsChooser(const std::string &pathToDir, choose_func_t chooseFunc) : _pathToDir(pathToDir),
                                                                                 _chooseFunc(chooseFunc) {}

        void choose();

    private:
        std::string _pathToDir;
        choose_func_t _chooseFunc;
    };

    std::vector<int> correlationColumnsChoose(const std::string &bucketPath);


} // namespace qtr