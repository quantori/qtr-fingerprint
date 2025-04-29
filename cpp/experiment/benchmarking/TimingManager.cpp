#include "TimingManager.h"

#include <glog/logging.h>

TimingManager::~TimingManager() {
    finish();
}

TimingManager::TimingManager(double timeLimitSeconds)
        : timeLimit(timeLimitSeconds),
          stopFlag(false),
          experimentFinished(false) {}

void TimingManager::start() {
    std::lock_guard<std::mutex> lock(mutex);
    stopFlag = false;
    experimentFinished.store(false, std::memory_order_release);
    startTime = std::chrono::high_resolution_clock::now();
    timeoutThread = std::thread(&TimingManager::timeoutWatcher, this);
}

double TimingManager::finish() {
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (experimentFinished.load(std::memory_order_acquire)) {
            return lastDuration;  // Already finished
        }
        auto now = std::chrono::high_resolution_clock::now();
        lastDuration = std::chrono::duration<double>(now - startTime).count();
        experimentFinished.store(true, std::memory_order_release);
    }

    if (timeoutThread.joinable()) {
        timeoutThread.join();
    }

    return lastDuration;
}

bool &TimingManager::getStopFlag() {
    return stopFlag;
}

void TimingManager::timeoutWatcher() {
    auto sleepDuration = std::chrono::duration<double>(
            std::clamp(timeLimit / 100.0, 0.001, 0.01)
    );

    while (!experimentFinished.load(std::memory_order_relaxed)) {
        std::this_thread::sleep_for(sleepDuration);
        if (checkTimeout()) {
            LOG(INFO) << "[TimingManager] Timeout triggered. Exiting watcher thread.";
            break;
        }
    }
}

bool TimingManager::checkTimeout() {
    std::lock_guard<std::mutex> lock(mutex);
    if (stopFlag) {
        return true;
    }

    auto now = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(now - startTime).count();
    if (elapsed > timeLimit) {
        stopFlag = true;
        LOG(INFO) << "[TimingManager] Timeout after " << elapsed << " seconds.";
        return true;
    }
    return false;
}

