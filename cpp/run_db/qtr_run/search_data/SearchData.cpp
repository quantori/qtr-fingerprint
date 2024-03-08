#include "SearchData.h"

namespace qtr {
    SearchData::SearchData(size_t ansCount, size_t threadCount, double timeLimit, bool verificationStage) :
            ansCount(ansCount), threadsCount(threadCount), timeLimit(timeLimit), verificationStage(verificationStage) {}
} // qtr
