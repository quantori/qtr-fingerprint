#pragma once

#include <vector>
#include <iostream>
#include <algorithm>

class ExperimentStat {
public:
    struct Entity {
        double duration;
        int resultsFound;
    };

    std::vector<Entity> _entities;

    const std::vector<Entity>& entities() const;

    ExperimentStat() = default;

    void add(const Entity &entity);

    [[nodiscard]] std::vector<double> Percentages(std::initializer_list<double> percetages) const;
};

std::ostream& operator << (std::ostream& out, const ExperimentStat& stat);
