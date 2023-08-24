#include "PreprocessingArgs.h"

using namespace std;

ABSL_FLAG(string, preprocessingType, "",
          "Source files type. "
          "Possible types: "
                  FLAG_NAME(SDF) ", "
                  FLAG_NAME(CSV));

ABSL_FLAG(string, sourceDir, "", "Directory with source files");

ABSL_FLAG(string, destDir, "", "Destination directory where preprocessed files should be stored");

ABSL_FLAG(string, targetFilesType, "",
          "Target files type (for " FLAG_NAME(SDF) " preprocessing type only). "
                                                   "Possible types: "
                  FLAG_NAME(RawBucket) ", "
                  FLAG_NAME(Tables));
