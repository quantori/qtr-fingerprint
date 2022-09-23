#include "gtest/gtest.h"

#include "../utils/TmpDirFixture.h"

#include "Fingerprint.h"
#include "BallTree.h"
#include "split_bit_selection/MaxDispersionBitSelector.h"
#include "io/BufferedWriter.h"
#include "io/BufferedReader.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"
#include "fingerprint_table_io/FingerprintTableReader.h"

class HighLevelBallTreeTests : public TmpDirFixture {
public:

    using DataTable = std::vector<std::pair<uint64_t, qtr::IndigoFingerprint>>;
    using QueryTable = std::vector<qtr::IndigoFingerprint>;

    size_t _treeDepth = 0;
    size_t _parallelizeDepth = 0;
    size_t _drivesNumber = 0;
    size_t _startSearchDepth = 0;
    bool _propertiesInitialized = false;

    void initProperties(size_t treeDepth = 0, size_t parallelizeDepth = 0, size_t drivesNumber = 0,
                        size_t startSearchDepth = 0) {
        if (_propertiesInitialized)
            GTEST_FAIL() << "initialize properties more the once per test";
        _treeDepth = treeDepth;
        if (parallelizeDepth >= treeDepth || startSearchDepth > treeDepth) {
            GTEST_FAIL() << "Wrong test arguments";
        }
        _parallelizeDepth = parallelizeDepth;
        _drivesNumber = drivesNumber;
        _startSearchDepth = startSearchDepth;
        _propertiesInitialized = true;
    }

    static QueryTable transformDataToQueries(const DataTable &data) {
        QueryTable result;
        for (const auto &[id, fp]: data) {
            result.emplace_back(fp);
        }
        return result;
    }

    static DataTable getAll5BitsMasksData() {
        DataTable result;
        for (size_t i = 0; i < 32; i++) {
            qtr::IndigoFingerprint fp;
            fp.reset();
            fp[0] = i & 1;
            fp[1] = i & 2;
            fp[2] = i & 4;
            fp[3] = i & 8;
            fp[4] = i & 16;
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

    static QueryTable getAll5BitsMasksQueries() {
        return transformDataToQueries(getAll5BitsMasksData());
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

    static QueryTable getRowOfSubMasksQueries() {
        return transformDataToQueries(getRowOfSubMasksData());
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

    static QueryTable getNoSubMasksQueries() {
        return transformDataToQueries(getNoSubMasksData());
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

    static QueryTable getHalfZerosAndHalfOnesQueries() {
        QueryTable result;
        qtr::IndigoFingerprint fingerprint;
        fingerprint.reset();
        result.emplace_back(fingerprint);
        fingerprint.set();
        result.emplace_back(fingerprint);
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

    static QueryTable getAllQueries() {
        return transformDataToQueries(getAllData());
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
        for (size_t i = 1; i <= _drivesNumber; i++) {
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

    void buildBallTree(const DataTable &data) {
        LOG(INFO) << "Start data preparation";
        prepareTreeDirs(data);
        LOG(INFO) << "Finish data preparation";
        LOG(INFO) << "Start building ball tree";
        qtr::BallTree ballTree(_treeDepth, _parallelizeDepth, getTreeDirs(), qtr::MaxDispersionBitSelector());
        LOG(INFO) << "Finish building ball tree";
        LOG(INFO) << "Start dumping ball tree to " << getTreePath();
        qtr::BufferedWriter treeWriter(getTreePath());
        ballTree.dumpNodes(treeWriter);
        LOG(INFO) << "Finish dumping ball tree to " << getTreePath();
    }

    void buildBallTreeAndCheck(const DataTable &data) {
        if (!_propertiesInitialized) {
            GTEST_FAIL() << "Properties wasn't initialized";
        }
        std::map<uint64_t, qtr::IndigoFingerprint> dataMap(data.begin(), data.end());
        std::set<uint64_t> dataUsed;
        buildBallTree(data);
        for (auto &dataDir: getTreeDirs()) {
            for (auto &leafDir: qtr::findFiles(dataDir, "")) {
                auto leafDataFile = leafDir / ("data" + qtr::fingerprintTableExtension);
                EXPECT_TRUE(std::filesystem::exists(leafDataFile));
                for (const auto &[id, fingerprint]: qtr::FingerprintTableReader(leafDataFile)) {
                    EXPECT_FALSE(dataUsed.count(id));
                    dataUsed.insert(id);
                    EXPECT_TRUE(dataMap.count(id));
                    const auto &[expectedId, expectedFingerprint] = *dataMap.find(id);
                    for (size_t i = 0; i < qtr::IndigoFingerprint::size(); i++) {
                        EXPECT_EQ(expectedFingerprint[i], fingerprint[i]);
                    }
                }
            }
        }
        EXPECT_EQ(dataUsed.size(), dataMap.size());
    }

    void runTest(const DataTable &data, const QueryTable &queries) {
        if (!_propertiesInitialized) {
            GTEST_FAIL() << "Properties wasn't initialized";
        }
        buildBallTreeAndCheck(data);
        qtr::BufferedReader treeReader(getTreePath());
        qtr::BallTree ballTree(treeReader, getTreeDirs());
        for (const auto &query: queries) {
            auto expectedAnswer = getAnswers(data, query);
            auto actualAnswer = ballTree.search(query, -1, _startSearchDepth);
            EXPECT_TRUE(std::is_permutation(expectedAnswer.begin(), expectedAnswer.end(),
                                            actualAnswer.begin(), actualAnswer.end()));
        }
    }
};

TEST_F(HighLevelBallTreeTests, buildBallTree) {
    initProperties(4, 2, 4, 2);
    buildBallTreeAndCheck(getAllData());
}

TEST_F(HighLevelBallTreeTests, buildANdRunBallTree5BitsMasks) {
    initProperties(6, 3, 4, 3);
    runTest(getAllData(), getAll5BitsMasksQueries());
}

TEST_F(HighLevelBallTreeTests, buildANdRunBallTreeRowOfSubMasks) {
    initProperties(6, 3, 4, 3);
    runTest(getAllData(), getRowOfSubMasksQueries());
}

TEST_F(HighLevelBallTreeTests, buildANdRunBallTreeNoSubMasks) {
    initProperties(6, 3, 4, 3);
    runTest(getAllData(), getNoSubMasksQueries());
}

TEST_F(HighLevelBallTreeTests, buildANdRunBallTreeAllZerosAndAllOnesQueries) {
    initProperties(6, 3, 4, 3);
    runTest(getAllData(), getHalfZerosAndHalfOnesQueries());
}