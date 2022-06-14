#pragma once

#include <cstdint>
#include <map>
#include <vector>

namespace qtr {

template<class Table>
class TableView {
public:
    using IndexType = uint32_t;

    explicit TableView(const Table *table = nullptr)
        : _table(table)
    {
        if (!table)
            return;

        _indices.resize(table->size());
        for(IndexType i = 0; i < _indices.size(); i++)
            _indices[i] = i;
    }


    template<class Filter>
    TableView<Table> filter(const Filter &filter) const {
        
        TableView<Table> result;
        result._table = _table;

        for(IndexType index : _indices) {
            if (filter(_table->at(index)))
                result._indices.push_back(index);
        }

        return result;
    }

    template<class Splitter>
    std::map<std::size_t, TableView<Table>> split(const Splitter &splitter) const {
        
        std::map<std::size_t, TableView<Table>> result;

        for(IndexType index : _indices) {
            std::size_t part = splitter(_table->at(index));
            result[part]._indices.push_back(index);
        }

        for(std::pair<const std::size_t, TableView<Table>> &pair: result) {
            pair.second._table = _table;
        }

        return result;
    }

    std::vector<IndexType>::const_iterator begin() const { return _indices.cbegin(); }
    std::vector<IndexType>::const_iterator end() const { return _indices.cend(); }


private:
    std::vector<IndexType> _indices;
    const Table *_table;
};

} // namespace qtr