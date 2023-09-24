#include "PreprocessorFactory.h"
#include "Utils.h"

using namespace qtr;
using namespace std;

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "preprocessing.info", true);
    PreprocessingArgs args(argc, argv);
    auto preprocessor = PreprocessorFactory::create(args.preprocessingType());
    preprocessor->run(args);
    return 0;
}