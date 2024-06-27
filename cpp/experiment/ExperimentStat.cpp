#include "ExperimentStat.h"

#include <iomanip>

const std::vector<ExperimentStat::Entity> &ExperimentStat::entities() const {
    return _entities;
}

void ExperimentStat::add(const ExperimentStat::Entity &entity) {
    _entities.push_back(entity);
}

std::vector<double> ExperimentStat::Percentages(std::initializer_list<double> percetages) const {
    auto entities = _entities;
    std::sort(entities.begin(), entities.end(), [](Entity x, Entity y) {
        return x.duration < y.duration;
    });
    std::vector<double> res;
    for (auto &p: percetages) {
        res.push_back(entities[size_t(p * (double) entities.size())].duration);
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
