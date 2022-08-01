#include "ColumnsChooser.h"
#include "Utils.h"
#include <fstream>
#include <Fingerprint.h>
#include <assert.h>
#include <math.h>
#include <numeric>
#include <future>

namespace qtr {

    /**
    * Call func for each record in data file
    * @tparam Functor call functor on IndigoFingerprint and string that represent smiles
    * @param filename
    * @param readSmiles if false than smiles equals to an empty string
    * @param func
    */
    template<typename Functor>
    void forEachLine(const std::string &filename, bool readSmiles, Functor &&func) {
        std::ifstream input(filename);
        if (input.bad()) {
            std::cerr << "smth with " << filename << " file\n";
            return;
        }
        uint64_t cntRecords;
        input.read((char *) (&cntRecords), sizeof(cntRecords));
        for (uint64_t record = 0; record < cntRecords; ++record) {
            IndigoFingerprint currFP;
            // Read fingerprint
            currFP.readFrom(input);
            // Read smiles
            if (readSmiles) { // TODO buffer reading?
                std::string smiles = "";
                char c;
                while ((c = input.get()) != '\n')
                    smiles += c;
                func(currFP, smiles);
            } else {
                while (input.get() != '\n') {}
                func(currFP, "");
            }
        }
    }

    void saveColumns(const std::string &bucketFileName, const std::vector<int> &columns) {
        std::ofstream out(bucketFileName + "OrderColumns.txt");
        for (auto col: columns)
            out << col << ' ';
    }

    void chooseAndSave(const std::string &bucketFileName, choose_func_t chooseFunc) {
        auto columns = chooseFunc(bucketFileName);
        saveColumns(bucketFileName, columns);
    }


    void ColumnsChooser::choose() {
//        #pragma omp parallel for
//        for (auto bucketFilePath : findFiles(_pathToDir, "")) {
//            std::cerr << "begin: " << bucketFilePath << std::endl;
//            chooseAndSave(bucketFilePath, _chooseFunc);
//            std::cerr << "  end: " << bucketFilePath << std::endl;
//        }
        using future_t = decltype(std::async(std::launch::async, chooseAndSave, "", _chooseFunc));
        std::vector<future_t> tasks;
        int started = 0;
        for (auto bucketFilePath: findFiles(_pathToDir, "")) {
            tasks.emplace_back(std::async(std::launch::async, chooseAndSave, bucketFilePath, _chooseFunc));
            std::cerr << "start: " << ++started << std::endl;
        }
        int completed = 0;
        for (auto &task: tasks) {
            task.get();
            std::cerr << "complete: " << ++completed << std::endl;
        }
    }

    double findCorrelation(const std::vector<bool> &x, const std::vector<bool> &y, uint32_t xSum, uint32_t ySum) {
        assert(x.size() == y.size());
        size_t c[2][2] = {0, 0, 0, 0};
        for (size_t i = 0; i < x.size(); i++) {
            c[x[i]][y[i]]++;
        }
        double xMean = xSum / x.size();
        double yMean = ySum / y.size();
        double numerator = 0;
        for (size_t i = 0; i <= 1; i++) {
            for (size_t j = 0; j <= 1; j++) {
                numerator += c[i][j] * (i - xMean) * (j - yMean);
            }
        }
        double denominator = (xSum * (1 - xMean) * (1 - xMean) + (1 - xSum) * xMean * xMean) *
                             (ySum * (1 - yMean) * (1 - yMean) + (1 - yMean) * yMean * yMean);
        return numerator / sqrt(denominator);
    }

    std::vector<int> correlationColumnsChoose(const std::string &bucketPath) {
        const size_t fingerprintSize = fromBytesToBits(IndigoFingerprint::sizeInBytes);
        std::vector<std::vector<bool>> columns(fingerprintSize);
        std::vector<uint32_t> columnSum(fingerprintSize);
        forEachLine(bucketPath, false, [&columns, &columnSum](const IndigoFingerprint &fp, const std::string &_) {
            for (size_t i = 0; i < fingerprintSize; i++) {
                columnSum[i] += fp[i];
                columns[i].emplace_back(fp[i]);
            }
        });
        assert(!columns[0].empty());
        std::vector<double> maxCorrelation(fingerprintSize);
        for (size_t i = 0; i < fingerprintSize; i++) {
            for (size_t j = i + 1; j < fingerprintSize; j++) {
                double corr = abs(findCorrelation(columns[i], columns[j], columnSum[i], columnSum[j]));
                maxCorrelation[i] = std::max(maxCorrelation[i], corr);
                maxCorrelation[j] = std::max(maxCorrelation[j], corr);
            }
        }
        std::vector<int> columnsOrder(fingerprintSize);
        std::iota(columnsOrder.begin(), columnsOrder.end(), 0);
        std::sort(columnsOrder.begin(), columnsOrder.end(), [&maxCorrelation](int a, int b) {
            return maxCorrelation[a] < maxCorrelation[b];
        });
        return columnsOrder;
    }

} // namespace qtr