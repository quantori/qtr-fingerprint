#include "SmilesDirParser.h"

#include <algorithm>
#include <execution>
#include <vector>
#include <string>

#include <glog/logging.h>


#include "SmilesCSVParser.h"


SmilesDirParser::SmilesDirParser(std::filesystem::path dirPath) : _dirPath(std::move(dirPath)) {
}

SmilesStorage SmilesDirParser::parse() const {
    auto dirIt = std::filesystem::directory_iterator(_dirPath);
    std::mutex mutex;
    std::vector<std::string> result;
    std::for_each(std::execution::par_unseq, std::filesystem::begin(dirIt), std::filesystem::end(dirIt),
                  [&](const std::filesystem::directory_entry &entry) {
                      auto &filename = entry.path();
                      if (!std::filesystem::is_regular_file(filename)) {
                          LOG(WARNING) << "Skipped dataset file: " << filename;
                          return;
                      }
                      auto parser = SmilesCSVParser(filename);
                      parser.parse(result, mutex);
                  });
    return SmilesStorage(std::move(result));
}


