#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <csignal>
#include <memory>
#include <filesystem>

#include "ServiceManager.h"
#include "WebService.h"
#include "io/SmilesDirParser.h"
#include "search/engines/SearchEngineType.h"
#include "search/engines/RDKitSearchEngine.h"
#include "search/engines/IndigoSearchEngine.h"
#include "search/engines/BallTreeSearchEngine.h"
#include "frameworks/RDKitFramework.h"
#include "frameworks/IndigoFramework.h"
#include "utils/Config.h"

// Command line flags
ABSL_FLAG(std::string, engine, "RDKit", "Search engine type (RDKit, Indigo, BallTreeRDKit, BallTreeIndigo)");
ABSL_FLAG(std::string, dataset, "", "Path to dataset directory containing CSV files");
ABSL_FLAG(int, port, 8080, "Port number for the web service");
ABSL_FLAG(int, depth, -1, "Ball tree depth (for BallTree engines, -1 for auto)");
ABSL_FLAG(double, fpRatio, 1.0, "Fingerprint ratio (for BallTree engines)");

// Signal handling flag
volatile std::sig_atomic_t g_shutdown = 0;

void signalHandler(int signal) {
    LOG(INFO) << "Received signal " << signal << ", shutting down gracefully...";
    g_shutdown = 1;
}

SearchEngineType parseEngineType(const std::string& engineStr) {
    if (engineStr == "RDKit") return SearchEngineType::RDKit;
    if (engineStr == "Indigo") return SearchEngineType::Indigo;
    if (engineStr == "BallTreeRDKit") return SearchEngineType::BallTreeRDKit;
    if (engineStr == "BallTreeIndigo") return SearchEngineType::BallTreeIndigo;
    
    throw std::invalid_argument("Unknown engine type: " + engineStr);
}

template<typename SearchEngineT>
void runService(Config config, const std::string& datasetPath, int port) {
    // Create service manager
    auto serviceManager = std::make_shared<ServiceManager<SearchEngineT>>(config);
    
    // Load dataset
    LOG(INFO) << "Loading dataset from: " << datasetPath;
    SmilesDirParser parser(datasetPath);
    auto dataset = parser.parse();
    LOG(INFO) << "Loaded " << dataset.size() << " molecules";
    
    // Initialize search engine
    LOG(INFO) << "Initializing search engine...";
    serviceManager->initializeEngine(std::move(dataset));
    LOG(INFO) << "Search engine ready";
    
    // Create and configure web service
    auto webService = std::make_unique<WebService<SearchEngineT>>(serviceManager, port);
    
    // Setup signal handlers for graceful shutdown
    std::signal(SIGINT, signalHandler);
    std::signal(SIGTERM, signalHandler);
    
    // Start the web service (blocking call)
    LOG(INFO) << "Web service starting on port " << port;
    webService->run();
    
    LOG(INFO) << "Service shutdown complete";
}

int main(int argc, char* argv[]) {
    // Initialize logging
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "service.info");
    FLAGS_alsologtostderr = true;
    
    try {
        // Parse command line arguments
        absl::ParseCommandLine(argc, argv);
        
        auto engineStr = absl::GetFlag(FLAGS_engine);
        auto datasetPath = absl::GetFlag(FLAGS_dataset);
        auto port = absl::GetFlag(FLAGS_port);
        auto depth = absl::GetFlag(FLAGS_depth);
        auto fpRatio = absl::GetFlag(FLAGS_fpRatio);
        
        // Validate required parameters
        if (datasetPath.empty()) {
            LOG(ERROR) << "Dataset path is required. Use --dataset flag.";
            return 1;
        }
        
        if (!std::filesystem::exists(datasetPath)) {
            LOG(ERROR) << "Dataset path does not exist: " << datasetPath;
            return 1;
        }
        
        if (!std::filesystem::is_directory(datasetPath)) {
            LOG(ERROR) << "Dataset path is not a directory: " << datasetPath;
            return 1;
        }
        
        // Parse engine type
        SearchEngineType engineType;
        try {
            engineType = parseEngineType(engineStr);
        } catch (const std::exception& e) {
            LOG(ERROR) << e.what();
            LOG(ERROR) << "Available engine types: RDKit, Indigo, BallTreeRDKit, BallTreeIndigo";
            return 1;
        }
        
        LOG(INFO) << "Starting QTR Fingerprint Search Service";
        LOG(INFO) << "Engine: " << engineStr;
        LOG(INFO) << "Dataset: " << datasetPath;
        LOG(INFO) << "Port: " << port;
        
        // Setup configuration
        Config config;
        if (depth >= 0) {
            config.set("depth", std::to_string(depth));
        }
        config.set("fpRatio", std::to_string(fpRatio));
        
        // Run service with appropriate template instantiation
        switch (engineType) {
            case SearchEngineType::RDKit:
                runService<RDKitSearchEngine>(config, datasetPath, port);
                break;
            case SearchEngineType::Indigo:
                runService<IndigoSearchEngine>(config, datasetPath, port);
                break;
            case SearchEngineType::BallTreeRDKit:
                runService<BallTreeSearchEngine<RDKitFramework>>(config, datasetPath, port);
                break;
            case SearchEngineType::BallTreeIndigo:
                runService<BallTreeSearchEngine<IndigoFramework>>(config, datasetPath, port);
                break;
            default:
                LOG(ERROR) << "Unsupported search engine type";
                return 1;
        }
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Service failed: " << e.what();
        return 1;
    } catch (...) {
        LOG(ERROR) << "Unknown error occurred";
        return 1;
    }
    
    return 0;
} 