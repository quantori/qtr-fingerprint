#include "PreprocessingArgs.h"

using namespace std;

ABSL_FLAG(std::string, preprocessingType, "",
          "Source files type. "
          "Possible types: "
                  FLAG_NAME(SDF) ", "
                  FLAG_NAME(CSV));

ABSL_FLAG(std::string, sourceDir, "", "Directory with source files");

ABSL_FLAG(std::string, destDir, "", "Destination directory where preprocessed files should be stored");

ABSL_FLAG(std::string, targetFilesType, "",
          "Target files type (for " FLAG_NAME(SDF) " preprocessing type only). "
                                                   "Possible types: "
                  FLAG_NAME(RawBucket) ", "
                  FLAG_NAME(Tables));
