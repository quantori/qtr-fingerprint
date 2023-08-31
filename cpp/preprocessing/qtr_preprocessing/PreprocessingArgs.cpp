#include "PreprocessingArgs.h"

using namespace std;

ABSL_FLAG(string, preprocessingType, "",
          "Source files type. "
          "Possible types: "
                  FLAG_NAME(SDF) ", "
                  FLAG_NAME(CSV));

ABSL_FLAG(string, preprocessDir, "", "Directory with source files");

ABSL_FLAG(string, destDir, "", "Destination directory where preprocessed files should be stored");

ABSL_FLAG(string, targetFilesType, "",
          "Target files type (for " FLAG_NAME(SDF) " preprocessing type only). "
                                                   "Possible types: "
                  FLAG_NAME(RawBucket) ", "
                  FLAG_NAME(Tables));

ABSL_FLAG(bool, properties, true,
          "Is true if properties should be concerned. False otherwise (for " FLAG_NAME(
                  CSV) " preprocessing type only)");