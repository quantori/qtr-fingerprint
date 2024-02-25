#include <future>
#include "SDFPreprocessor.h"

#include "Utils.h"
#include "Fingerprint.h"
#include "raw_bucket_io/RawBucketWriter.h"
#include "string_table_io/StringTableWriter.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"

#include "IndigoWriteBuffer.h"
#include "IndigoIterator.h"
#include "Profiling.h"
#include "id_to_string_io/IdToStringWriter.h"

using namespace qtr;
using namespace std;

namespace {
    void processSDF(const filesystem::path &sdFilePath,
                    const function<void(indigo_cpp::IndigoMoleculeSPtr)> &processMol) {
        LOG(INFO) << "Start processing " << sdFilePath;
        auto indigoSessionPtr = indigo_cpp::IndigoSession::create();
        uint64_t skipped = 0;
        uint64_t processed = 0;
        for (auto &mol: indigoSessionPtr->iterateSDFile(sdFilePath)) {
            try {
                processMol(mol);
                processed++;
            }
            catch (const exception &e) {
                LOG(ERROR) << "Fail to parse molecule from " << sdFilePath << " -- " << e.what();
                skipped++;
            }
            if ((processed + skipped) % 100'000 == 0)
                LOG(INFO) << (processed + skipped) << " molecules was processed from " << sdFilePath;
        }
        LOG(INFO) << "Finish processing " << sdFilePath << " : skipped -- " << skipped << ", processed -- "
                  << processed;
    }

    void parseSDF(const filesystem::path &sdFilePath, const PreprocessingArgs &args) {
        string fileName = sdFilePath.stem();
        while (count(fileName.begin(), fileName.end(), '.') != 0)
            fileName.pop_back();
        filesystem::path smilesTablePath =
                args.smilesTables() / (fileName + stringTableExtension);
        filesystem::path fingerprintTablePath =
                args.fingerprintTables() / (fileName + fingerprintTableExtension);
        filesystem::path idToStringTablePath =
                args.idToStringTables() / (fileName + ".csv");


        filesystem::create_directory(smilesTablePath.parent_path());
        filesystem::create_directory(fingerprintTablePath.parent_path());
        filesystem::create_directory(idToStringTablePath.parent_path());

        StringTableWriter smilesTableWriter(smilesTablePath);
        FingerprintTableWriter fingerprintTableWriter(fingerprintTablePath);
        IdToStringWriter idToStringWriter(idToStringTablePath);

        static atomic_uint64_t entryId = 0;

        processSDF(sdFilePath,
                   [&smilesTableWriter, &fingerprintTableWriter, &idToStringWriter](
                           const indigo_cpp::IndigoMoleculeSPtr &mol) {
                       mol->aromatize();
                       string smiles = mol->smiles();
                       Fingerprint fingerprint = indigoFingerprintFromSmiles(smiles);
                       uint64_t cid = stoull(mol->name());

                       smilesTableWriter << make_pair(cid, smiles);
                       fingerprintTableWriter << make_pair(cid, fingerprint);
                       idToStringWriter << make_pair(entryId.fetch_add(1), to_string(cid));
                   });
    }

}

void SDFPreprocessor::run(const PreprocessingArgs &args) {
    ProfileScope("SDF preprocessing");
    vector<filesystem::path> sdFilePaths = findFiles(args.sourceDir(), ".gz");
    vector<future<void>> tasks;
    for (auto &sdFilePath: sdFilePaths) {
        tasks.emplace_back(async(launch::async, parseSDF, cref(sdFilePath), cref(args)));
    }
    for (auto &task: tasks) {
        task.get();
    }
}
