#pragma once

#include <future>
#include <vector>
#include <memory>

#include "Fingerprint.h"
#include "answer_filtering/AnswerFilter.h"
#include "answer_filtering/AlwaysTrueFilter.h"

namespace qtr {

    class BallTreeQueryData {

    public:
        explicit BallTreeQueryData(size_t stopAnswersCount, const IndigoFingerprint &query = IndigoFingerprint(),
                                   std::unique_ptr<AnswerFilter> &&filter = std::make_unique<AlwaysTrueFilter>());

        void addAnswers(const std::vector<CIDType> &answers);

        void filterAndAddAnswers(const std::vector<CIDType> &answers, AnswerFilter &filterObject);

        [[nodiscard]] bool isFinished() const;

        [[nodiscard]] bool checkFoundEnoughAnswers() const;

        void stopProcess();

        void waitAllTasks();

        std::pair<bool, std::vector<CIDType>> getAnswers(size_t beginPos, size_t endPos);

        void addTask(std::future<void> &&task);

        void tagFinishTask();

        void tagStartTask();

        [[nodiscard]] std::unique_ptr<AnswerFilter> getFilterObject() const;

        [[nodiscard]] size_t getCurrentAnswersCount() const;

        [[nodiscard]] const IndigoFingerprint &getQueryFingerprint() const;

        [[nodiscard]] bool checkShouldStop() const;

    private:
        IndigoFingerprint _queryFingerprint;
        std::vector<CIDType> _result;
        std::mutex _resultLock;
        size_t _stopAnswersNumber;
        std::atomic_size_t _startedTasksCount;
        std::atomic_size_t _finishedTasksCount;
        std::unique_ptr<AnswerFilter> _filter;
        std::atomic_bool _shouldStopProcess;
        std::vector<std::future<void>> _tasks;
    };

} // namespace qtr
