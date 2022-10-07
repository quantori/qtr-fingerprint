#include "indigo.h"

#include "IndigoMolecule.h"
#include "IndigoSession.h"
#include "IndigoWriteBuffer.h"
#include "IndigoSDFileIterator.h"

#include <glog/logging.h>

#include <chrono>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <future>
#include <functional>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include "Utils.h"
#include "Fingerprint.h"
#include "raw_bucket_io/RawBucketWriter.h"
#include "smiles_table_io/SmilesTableWriter.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"


ABSL_FLAG(std::string, source_dir_path, "",
          "Path to dir with sdf files");

ABSL_FLAG(std::string, dest_dir_path, "",
          "Path to directory where parsed data should be stored");

ABSL_FLAG(std::string, dest_type, "",
          R"(How parsed data should be stored: "rb", "tables")");

enum DestType {
    RB,
    TABLES
};

struct Args {
    std::filesystem::path sourceDirPath;
    std::filesystem::path destDirPath;
    DestType destType;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        sourceDirPath = absl::GetFlag(FLAGS_source_dir_path);
        qtr::emptyArgument(sourceDirPath, "Please specify source_dir_path option");

        destDirPath = absl::GetFlag(FLAGS_dest_dir_path);
        qtr::emptyArgument(destDirPath, "Please specify dest_dir_path option");

        std::string parseModeStr = absl::GetFlag(FLAGS_dest_type);
        if (parseModeStr == "rb") {
            destType = RB;
        } else if (parseModeStr == "tables") {
            destType = TABLES;
        } else {
            LOG(ERROR) << R"(Please specify parseModeStr option with value "rb" or "tables")";
            exit(-1);
        }
    }
};

void processSDF(const std::filesystem::path &sdFilePath,
                const std::function<void(indigo_cpp::IndigoMoleculeSPtr)> &processMol) {
    LOG(INFO) << "Start processing " << sdFilePath;
    auto indigoSessionPtr = indigo_cpp::IndigoSession::create();
    uint64_t skipped = 0;
    uint64_t processed = 0;
    for (auto &mol: indigoSessionPtr->iterateSDFile(sdFilePath)) {
        try {
            processMol(mol);
            processed++;
        }
        catch (const std::exception &e) {
            LOG(ERROR) << "Fail to parse molecule from " << sdFilePath << " -- " << e.what();
            skipped++;
        }
        if ((processed + skipped) % 100'000 == 0)
            LOG(INFO) << (processed + skipped) << " molecules was processed from " << sdFilePath;
    }
    LOG(INFO) << "Finish processing " << sdFilePath << " : skipped -- " << skipped << ", processed -- " << processed;
}

void sdfToRb(const std::filesystem::path &sdFilePath, const Args &args) {
    std::filesystem::path rbFilePath = args.destDirPath / (sdFilePath.stem().string() + qtr::rawBucketExtension);
    qtr::RawBucketWriter writer(rbFilePath);

    processSDF(sdFilePath, [&writer](const indigo_cpp::IndigoMoleculeSPtr &mol) {
        mol->aromatize();
        qtr::IndigoFingerprint fingerprint = qtr::IndigoFingerprintFromSmiles(mol->smiles());
        writer << std::make_pair(mol->smiles(), fingerprint);
    });
}

void sdfToTables(const std::filesystem::path &sdFilePath, const Args &args) {
    std::filesystem::path smilesTablePath = args.destDirPath / (sdFilePath.stem().string() + qtr::smilesTableExtension);
    std::filesystem::path fingerprintTablePath =
            args.destDirPath / (sdFilePath.stem().stem().string() + qtr::fingerprintTableExtension);

    qtr::SmilesTableWriter smilesTableWriter(smilesTablePath);
    qtr::FingerprintTableWriter fingerprintTableWriter(fingerprintTablePath);

    processSDF(sdFilePath, [&smilesTableWriter, &fingerprintTableWriter](const indigo_cpp::IndigoMoleculeSPtr &mol) {
        mol->aromatize();
        std::string smiles = mol->smiles();
        qtr::IndigoFingerprint fingerprint = qtr::IndigoFingerprintFromSmiles(smiles);
        uint64_t cid = std::stoull(mol->name());

        smilesTableWriter << std::make_pair(cid, smiles);
        fingerprintTableWriter << std::make_pair(cid, fingerprint);
    });
}

void parseSDF(const std::filesystem::path &sdFilePath, const Args &args) {
    if (args.destType == RB) {
        sdfToRb(sdFilePath, args);
    } else if (args.destType == TABLES) {
        sdfToTables(sdFilePath, args);
    }
}

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::LogToStderr();
    Args args(argc, argv);

    qtr::TimeTicker timeTicker;
    std::vector<std::filesystem::path> sdFilePaths = qtr::findFiles(args.sourceDirPath, ".sdf");
    std::vector<std::future<void>> tasks;
    for (auto &sdFilePath: sdFilePaths) {
        tasks.emplace_back(std::async(std::launch::async, parseSDF, std::cref(sdFilePath), std::cref(args)));
    }
    for (auto &task: tasks) {
        task.get();
    }
    timeTicker.tick("Elapsed time");
    return 0;
}