#pragma once

#include "Fingerprint.h"

namespace qtr {
    class AnswerFilter {
    public:
        virtual bool operator()(const IndigoFingerprint &fingerprint) = 0;

        virtual static AnswerFilter prepare() = 0;
    };
}