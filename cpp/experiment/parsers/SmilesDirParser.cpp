#include "SmilesDirParser.h"

#include "SmilesCSVParser.h"

#include <glog/logging.h>

#include <algorithm>
#include <execution>

SmilesDirParser::SmilesDirParser(std::filesystem::path dirPath) : _dirPath(std::move(dirPath)) {
}

void SmilesDirParser::parse(std::vector<std::string> &dest, std::mutex &m) const {
    auto dirIt = std::filesystem::directory_iterator(_dirPath);
    std::for_each(std::execution::par_unseq, std::filesystem::begin(dirIt), std::filesystem::end(dirIt),
                  [&](const std::filesystem::directory_entry &entry) {
                      auto &filename = entry.path();
                      if (!std::filesystem::is_regular_file(filename)) {
                          LOG(WARNING) << "Skipped dataset file: " << filename;
                          return;
                      }
                      auto parser = SmilesCSVParser(filename);
                      parser.parse(dest, m);
                  });
}


