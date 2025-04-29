#pragma once

#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>
#include <functional>
#include <iostream>

class TimingManager {
public:
    explicit TimingManager(double timeLimitSeconds);

    ~TimingManager();

    void start();

    double finish();

    [[nodiscard]] bool& getStopFlag();

private:
    double timeLimit;
    bool stopFlag;
    std::atomic<bool> experimentFinished;
    std::chrono::high_resolution_clock::time_point startTime;
    std::thread timeoutThread;
    std::mutex mutex;
    double lastDuration = 0.0;

    void timeoutWatcher();

    bool checkTimeout();
};
