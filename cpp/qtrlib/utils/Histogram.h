#pragma once

#include <cstdint>
#include <vector>

namespace qtr {

class Histogram {
public:
    using CounterType = uint32_t;
    
    Histogram(size_t binsCount) : _bins(binsCount, 0) {}
    void add(size_t bin, CounterType num) { _bins.at(bin) += num; }
    const std::vector<CounterType> &bins() const { return _bins; }

private:
    std::vector<CounterType> _bins;
};

} // namespace qtr