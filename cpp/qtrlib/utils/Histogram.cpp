#include "Histogram.h"

namespace qtr {

std::size_t Histogram::findClosestBin(CounterType value) const
{
    std::size_t result = std::size_t(-1);
    Histogram::CounterType deviation = value;

    for(std::size_t bin = 0; bin < _bins.size(); bin++) {
        
        CounterType num = _bins.at(bin);
        CounterType dev = (num > value ? num - value : value - num);
        
        if (dev < deviation) {
            deviation = dev;
            result = bin;
        }
    }

    return result;
}

} // namespace qtr