#include "PreprocessorFactory.h"
#include "Utils.h"
#include "PreprocessingArgs.h"

using namespace qtr;
using namespace std;

ABSL_FLAG(std::string, preprocessingType, "",
          "Source files type. "
          "Possible types: "
                  FLAG_NAME(SDF) ", "
                  FLAG_NAME(CSV));

ABSL_FLAG(std::string, sourceDir, "", "Directory with source files");

ABSL_FLAG(std::string, destDir, "", "Destination directory where preprocessed files should be stored");

ABSL_FLAG(std::string, targetFilesType, "",
          "Target files type. "
          "Possible types: "
                  FLAG_NAME(RawBucket) ", "
                  FLAG_NAME(Tables));


int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "preprocessing.info", true);
    PreprocessingArgs args(argc, argv);
    auto preprocessor = PreprocessorFactory::create(args.preprocessingType());
    preprocessor->run(args);
    return 0;
}