#include "gtest/gtest.h"

#include "../utils/TmpDirFixture.h"

#include "Fingerprint.h"
#include "BallTree.h"
#include "split_bit_selection/MaxDispersionBitSelector.h"
#include "io/BufferedWriter.h"
#include "io/BufferedReader.h"
#include "data_io/fingerprint_table_io/FingerprintTableWriter.h"

class HighLevelBallTreeTests : public TmpDirFixture {
public:

    using DataTable = std::vector<std::pair<uint64_t, qtr::IndigoFingerprint>>;
    using QueryTable = std::vector<qtr::IndigoFingerprint>;

    static DataTable getAll5BitsMasksData() {
        DataTable result;
        for (size_t i = 0; i < 32; i++) {
            qtr::IndigoFingerprint fp;
            fp.reset();
            fp[0] = i & 1;
            fp[1] = i & 2;
            fp[3] = i & 4;
            fp[4] = i & 8;
            fp[5] = i & 16;
            result.emplace_back(i, fp);
        }
        return result;
    }

    static DataTable getAll5BitMasksDataWithDuplicates() {
        auto result = HighLevelBallTreeTests::getAll5BitsMasksData();
        auto result2 = result;
        result.insert(result.end(), result2.begin(), result2.end());
        return result;
    }

    static DataTable getRowOfSubMasksData() {
        DataTable result;
        qtr::IndigoFingerprint fp;
        fp.reset();
        for (size_t i = 0; i < qtr::IndigoFingerprint::size(); i++) {
            fp[i] = true;
            result.emplace_back(i, fp);
        }
        return result;
    }

    static DataTable getNoSubMasksData() {
        DataTable result;
        qtr::IndigoFingerprint fp;
        for (size_t i = 0; i < qtr::IndigoFingerprint::size(); i++) {
            fp.reset();
            fp[i] = true;
            result.emplace_back(i, fp);
        }
        return result;
    }

    static DataTable getHalfZerosAndHalfOnesData() {
        DataTable result;
        qtr::IndigoFingerprint fp;
        fp.reset();
        for (size_t i = 0; i < 10; i++) {
            result.emplace_back(i, fp);
        }
        fp.set();
        for (size_t i = 10; i < 20; i++) {
            result.emplace_back(i, fp);
        }
        return result;
    }

    static DataTable getAllData() {
        DataTable result;
        auto all5BitsMasksWithDuplicates = getAll5BitMasksDataWithDuplicates();
        result.insert(result.end(), all5BitsMasksWithDuplicates.begin(), all5BitsMasksWithDuplicates.end());
        auto rowOfSubMasks = getRowOfSubMasksData();
        result.insert(result.end(), rowOfSubMasks.begin(), rowOfSubMasks.end());
        auto noSubMasks = getNoSubMasksData();
        result.insert(result.end(), noSubMasks.begin(), noSubMasks.end());
        auto halfZerosHalfOnes = getHalfZerosAndHalfOnesData();
        result.insert(result.end(), halfZerosHalfOnes.begin(), halfZerosHalfOnes.end());
        for (size_t i = 0; i < result.size(); i++) {
            result[i].first = i;
        }
        return result;
    }

    static QueryTable getQueries() {
        QueryTable result;
        for (const auto &[id, fp]: getAllData()) {
            result.emplace_back(fp);
        }
        return result;
    }

    static std::vector<uint64_t> getAnswers(const DataTable &data, const qtr::IndigoFingerprint &query) {
        std::vector<size_t> result;
        for (const auto &[id, fp]: data) {
            if (query <= fp) {
                result.emplace_back(id);
            }
        }
        return result;
    }

    std::filesystem::path getTreePath() const {
        return getTmpDir() / "tree";
    }

    std::vector<std::filesystem::path> getTreeDirs() const {
        std::vector<std::filesystem::path> result;
        for (size_t i = 1; i <= 4; i++) {
            auto dir = getTmpDir() / std::to_string(i);
            std::filesystem::create_directory(dir);
            result.emplace_back(dir);
        }
        return result;
    }

    void prepareTreeDirs(const DataTable &data) {
        auto treeDirs = getTreeDirs();
        for (size_t i = 0; i < treeDirs.size(); i++) {
            auto writerPath = treeDirs[i] / "0" / "_data.ft";
            LOG(INFO) << "Create dir: " << writerPath.parent_path();
            std::filesystem::create_directory(writerPath.parent_path());
            qtr::FingerprintTableWriter writer(writerPath);
            for (size_t j = i; j < data.size(); j += treeDirs.size()) {
                writer << data[j];
            }
        }
    }

    void buildBallTree(const DataTable &data, size_t depth, size_t parallelize_depth) {
        LOG(INFO) << "Start data preparation";
        prepareTreeDirs(data);
        LOG(INFO) << "Finish data preparation";
        LOG(INFO) << "Start building ball tree";
        qtr::BallTree ballTree(depth, parallelize_depth, getTreeDirs(), qtr::MaxDispersionBitSelector());
        LOG(INFO) << "Finish building ball tree";
        LOG(INFO) << "Start dumping ball tree to " << getTreePath();
        qtr::BufferedWriter<4096> treeWriter(getTreePath());
        ballTree.dumpNodes(treeWriter);
        LOG(INFO) << "Finish dumping ball tree to " << getTreePath();
    }

    void runTest(const DataTable &data, const QueryTable &queries, size_t depth, size_t parallelize_depth,
                 size_t startSearchDepth) {
        if (parallelize_depth >= depth || startSearchDepth > depth) {
            GTEST_FAIL() << "Wrong test arguments";
        }
        buildBallTree(data, depth, parallelize_depth);
        qtr::BufferedReader<4096> treeReader(getTreePath());
        qtr::BallTree ballTree(treeReader, getTreeDirs());
        for (const auto &query: queries) {
            auto expectedAnswer = getAnswers(data, query);
            auto actualAnswer = ballTree.search(query, -1, startSearchDepth);
            EXPECT_TRUE(std::is_permutation(expectedAnswer.begin(), expectedAnswer.end(),
                                            actualAnswer.begin(), actualAnswer.end()));
        }
    }
};

TEST_F(HighLevelBallTreeTests, AllDataTest) {
    runTest(getAll5BitsMasksData(), getQueries(), 4, 2, 2);
}