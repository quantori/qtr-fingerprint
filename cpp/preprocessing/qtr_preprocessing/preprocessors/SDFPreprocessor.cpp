#include <future>
#include "SDFPreprocessor.h"

#include "PreprocessingArgs.h"
#include "Utils.h"
#include "Fingerprint.h"
#include "raw_bucket_io/RawBucketWriter.h"
#include "string_table_io/StringTableWriter.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"

#include "IndigoMolecule.h"
#include "IndigoWriteBuffer.h"
#include "IndigoIterator.h"
#include "Profiling.h"


namespace qtr {
    namespace {
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
            LOG(INFO) << "Finish processing " << sdFilePath << " : skipped -- " << skipped << ", processed -- "
                      << processed;
        }

        void sdfToRb(const std::filesystem::path &sdFilePath, const PreprocessingArgs &args) {
            std::filesystem::path rbFilePath = args.destDir() / (sdFilePath.stem().string() + qtr::rawBucketExtension);
            qtr::RawBucketWriter writer(rbFilePath);

            processSDF(sdFilePath, [&writer](const indigo_cpp::IndigoMoleculeSPtr &mol) {
                mol->aromatize();
                qtr::IndigoFingerprint fingerprint = qtr::indigoFingerprintFromSmiles(mol->smiles());
                writer << std::make_pair(mol->smiles(), fingerprint);
            });
        }

        void sdfToTables(const std::filesystem::path &sdFilePath, const PreprocessingArgs &args) {
            std::filesystem::path smilesTablePath =
                    args.smilesTables() / (sdFilePath.stem().string() + qtr::stringTableExtension);
            std::filesystem::path fingerprintTablePath =
                    args.fingerprintTables() / (sdFilePath.stem().string() + qtr::fingerprintTableExtension);

            std::filesystem::create_directory(smilesTablePath.parent_path());
            std::filesystem::create_directory(fingerprintTablePath.parent_path());

            qtr::StringTableWriter smilesTableWriter(smilesTablePath);
            qtr::FingerprintTableWriter fingerprintTableWriter(fingerprintTablePath);

            processSDF(sdFilePath,
                       [&smilesTableWriter, &fingerprintTableWriter](const indigo_cpp::IndigoMoleculeSPtr &mol) {
                           mol->aromatize();
                           std::string smiles = mol->smiles();
                           qtr::IndigoFingerprint fingerprint = qtr::indigoFingerprintFromSmiles(smiles);
                           uint64_t cid = std::stoull(mol->name());

                           smilesTableWriter << std::make_pair(cid, smiles);
                           fingerprintTableWriter << std::make_pair(cid, fingerprint);
                       });
        }

        void parseSDF(const std::filesystem::path &sdFilePath, const PreprocessingArgs &args) {
            if (args.targetFilesType() == TargetType::RawBucket) {
                sdfToRb(sdFilePath, args);
            } else if (args.targetFilesType() == TargetType::Tables) {
                sdfToTables(sdFilePath, args);
            }
        }

    }

    void SDFPreprocessor::run(const PreprocessingArgs &args) {
        ProfileScope("SDF preprocessing");
        std::vector<std::filesystem::path> sdFilePaths = qtr::findFiles(args.preprocessDir(), ".sdf");
        std::vector<std::future<void>> tasks;
        for (auto &sdFilePath: sdFilePaths) {
            tasks.emplace_back(std::async(std::launch::async, parseSDF, std::cref(sdFilePath), std::cref(args)));
        }
        for (auto &task: tasks) {
            task.get();
        }
    }
} // qtr