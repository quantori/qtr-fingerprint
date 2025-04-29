#include "ExperimentStat.h"

#include <iomanip>
#include <cassert>

const std::vector<ExperimentStat::Entity> &ExperimentStat::entities() const {
    return _entities;
}

void ExperimentStat::add(const ExperimentStat::Entity &entity) {
    _entities.push_back(entity);
}

std::vector<double> ExperimentStat::quantiles(std::initializer_list<double> quantiles) const {
    auto entities = _entities;
    std::sort(entities.begin(), entities.end(), [](Entity x, Entity y) {
        return x.duration < y.duration;
    });
    std::vector<double> res;
    for (auto &p: quantiles) {
        auto idx = size_t(p * (double) (entities.size() - 1));
        assert(0 <= idx);
        assert(idx < entities.size());
        res.push_back(entities[idx].duration);
    }
    return res;
}

std::ostream &operator<<(std::ostream &out, const ExperimentStat &stat) {
    out << "Id,duration,results\n";
    const auto &entities = stat.entities();
    out.precision(5);
    out.setf(std::ios::fixed);
    for (size_t i = 0; i < entities.size(); i++) {
        out << i + 1
            << "," << entities[i].duration
            << "," << entities[i].resultsFound << "\n";
    }
    return out;
}
