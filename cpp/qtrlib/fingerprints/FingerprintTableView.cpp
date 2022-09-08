#include "FingerprintTableView.h"
#include "Histogram.h"

namespace qtr {

template<>
Histogram IndigoFingerprintTableView::histogram() const {
    Histogram histogram(CHAR_BIT*IndigoFingerprint::size);
    
    for(IndexType index : _indices) {
        const IndigoFingerprint &fp = _table->at(index);
        for(std::size_t bit = 0; bit < fp.size; bit++)
            histogram.add(bit, Histogram::CounterType(fp[bit]));
    }

    return histogram;
}

} // namespace qtr