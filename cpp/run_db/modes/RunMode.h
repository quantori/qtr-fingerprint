#pragma once

namespace qtr {
    class RunMode {
    public:
        virtual void run() = 0;

        virtual ~RunMode() = default;
    };
}