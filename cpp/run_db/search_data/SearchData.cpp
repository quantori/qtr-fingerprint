#include "SearchData.h"

namespace qtr {
    SearchData::SearchData(TimeTicker &timeTicker, size_t ansCount, size_t threadCount, double timeLimit) :
            timeTicker(timeTicker), ansCount(ansCount), threadsCount(threadCount), timeLimit(timeLimit) {}
} // qtr