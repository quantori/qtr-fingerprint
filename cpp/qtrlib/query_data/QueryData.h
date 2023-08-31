#pragma once

#include <atomic>
#include <vector>
#include <future>
#include <algorithm>

#include "BallTreeTypes.h"
#include "AnswerFilter.h"

namespace qtr {

    template<typename AnsType>
    class QueryData {
    public:
        inline static std::atomic_uint64_t timedOutCounter = 0;

        QueryData(size_t stopAnswersCount, double timeLimit, std::unique_ptr<AnswerFilter<AnsType>> &&filter);

        void addAnswers(const std::vector<AnsType> &answers);

        void filterAndAddAnswers(const std::vector<AnsType> &answers, AnswerFilter<AnsType> &filterObject);

        [[nodiscard]] bool checkFoundEnoughAnswers() const;

        [[nodiscard]] size_t getCurrentAnswersCount() const;

        void stopProcess();

        [[nodiscard]] bool checkTimeOut();

        [[nodiscard]] bool isFinished() const;

        void waitAllTasks();

        std::pair<bool, std::vector<AnsType>> getAnswers(size_t beginPos, size_t endPos);

        void addTask(std::future<void> &&task);

        void tagFinishTask();

        void tagStartTask();

        [[nodiscard]] std::unique_ptr<AnswerFilter<AnsType>> getFilterObject() const;

        [[nodiscard]] bool checkShouldStop() const;

        virtual ~QueryData();

    protected:
        size_t _stopAnswersNumber;
        double _timeLimit;
        std::mutex _resultLock;
        std::vector<AnsType> _result;
        std::atomic_bool _shouldStopProcess;
        std::atomic_bool _wasTimeOut;
        decltype(std::chrono::steady_clock::now()) _startTimePoint;
        std::atomic_size_t _startedTasksCount;
        std::atomic_size_t _finishedTasksCount;
        std::vector<std::future<void>> _tasks;
        std::unique_ptr<AnswerFilter<AnsType>> _filter;
    };

    template<typename AnsType>
    QueryData<AnsType>::~QueryData() {
        QueryData<AnsType>::timedOutCounter += _wasTimeOut;
    }

    template<typename AnsType>
    bool QueryData<AnsType>::checkShouldStop() const {
        return _shouldStopProcess;
    }

    template<typename AnsType>
    std::unique_ptr<AnswerFilter<AnsType>> QueryData<AnsType>::getFilterObject() const {
        return _filter->copy();
    }

    template<typename AnsType>
    void QueryData<AnsType>::tagStartTask() {
        _startedTasksCount++;
    }

    template<typename AnsType>
    void QueryData<AnsType>::tagFinishTask() {
        _finishedTasksCount++;
    }

    template<typename AnsType>
    void QueryData<AnsType>::addTask(std::future<void> &&task) {
        tagStartTask();
        _tasks.emplace_back(std::move(task));
    }

    template<typename AnsType>
    std::pair<bool, std::vector<AnsType>> QueryData<AnsType>::getAnswers(size_t beginPos, size_t endPos) {
        std::lock_guard<std::mutex> lock(_resultLock);
        beginPos = std::min(beginPos, _result.size());
        endPos = std::min(endPos, _result.size());
        std::vector<AnsType> result(_result.begin() + beginPos, _result.begin() + endPos);
        return {isFinished(), result};
    }

    template<typename AnsType>
    void QueryData<AnsType>::waitAllTasks() {
        for (auto &task: _tasks)
            task.wait();
    }

    template<typename AnsType>
    bool QueryData<AnsType>::isFinished() const {
        return _startedTasksCount == _finishedTasksCount;
    }

    template<typename AnsType>
    bool QueryData<AnsType>::checkTimeOut() {
        if (_wasTimeOut)
            return true;
        auto now = std::chrono::steady_clock::now();
        std::chrono::duration<double> duration = now - _startTimePoint;
        bool timeOut = duration.count() > _timeLimit;
        if (timeOut) {
            _wasTimeOut = true;
            _shouldStopProcess = true;
        }
        return timeOut;
    }

    template<typename AnsType>
    void
    QueryData<AnsType>::filterAndAddAnswers(const std::vector<AnsType> &answers, AnswerFilter<AnsType> &filterObject) {
        std::vector<AnsType> filteredAnswers;
        for (auto &ans: answers) {
            if (filterObject(ans)) {
                filteredAnswers.emplace_back(ans);
            }
            if (checkTimeOut())
                break;
        }
        if (filteredAnswers.empty())
            return;
        addAnswers(filteredAnswers);
    }

    template<typename AnsType>
    void QueryData<AnsType>::stopProcess() {
        _shouldStopProcess = true;
    }

    template<typename AnsType>
    size_t QueryData<AnsType>::getCurrentAnswersCount() const {
        return _result.size();
    }

    template<typename AnsType>
    bool QueryData<AnsType>::checkFoundEnoughAnswers() const {
        return getCurrentAnswersCount() >= _stopAnswersNumber;
    }

    template<typename AnswersNumber>
    void QueryData<AnswersNumber>::addAnswers(const std::vector<AnswersNumber> &answers) {
        std::lock_guard<std::mutex> lock(_resultLock);
        copy(answers.begin(), answers.end(), back_inserter(_result));
        if (checkFoundEnoughAnswers())
            stopProcess();
    }

    template<typename AnsType>
    QueryData<AnsType>::QueryData(size_t stopAnswersCount, double timeLimit,
                                  std::unique_ptr<AnswerFilter<AnsType>> &&filter) :
            _timeLimit(timeLimit), _stopAnswersNumber(stopAnswersCount), _filter(std::move(filter)),
            _startedTasksCount(0), _finishedTasksCount(0), _shouldStopProcess(false), _tasks(), _wasTimeOut(false),
            _startTimePoint(std::chrono::steady_clock::now()) {}

} // qtr
