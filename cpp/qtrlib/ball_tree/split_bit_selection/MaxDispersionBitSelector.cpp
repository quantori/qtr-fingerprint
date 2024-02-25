#include <cmath>
#include <future>

#include "MaxDispersionBitSelector.h"
#include "Fingerprint.h"

namespace qtr {

    namespace {
        size_t getBitValToMinimize(const ColumnsStatistic &statistic, size_t bit) {
            return (size_t) std::abs(2 * (int16_t) statistic.zeros(bit) - (int64_t) statistic.size());
        }
    }


    size_t MaxDispersionBitSelector::operator()(const ColumnsStatistic &statistic) const {
        size_t bestBit = 0;
        size_t bestBitVal = getBitValToMinimize(statistic, bestBit);
        for (size_t i = 1; i < statistic.columns(); i++) {
            size_t iVal = getBitValToMinimize(statistic, i);
            if (iVal < bestBitVal) {
                bestBit = i;
                bestBitVal = iVal;
            }
        }
        return bestBit;
    }

} // qtr