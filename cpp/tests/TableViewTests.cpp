#include "TableView.h"

#include "Common.h"

#include <gtest/gtest.h>

#include <vector>

namespace qtr {

class TestTable {
public:
    explicit TestTable(size_t n) {
        _data.resize(n);
        for(size_t i = 0; i < n; i++)
            _data.at(i) = int(i);
    }

    size_t size() const { return _data.size(); } 
    int at(size_t i) const { return _data.at(i); }

    const std::vector<int> &data() const { return _data; }

private:
    std::vector<int> _data;
};

} // namespace qtr

using namespace qtr;

TEST(TableView, SIMPLE_TABLE) {
    TestTable testTable(5);
    
    std::vector<TableView<TestTable>::IndexType> oddIndicesTest = {1, 3};
    std::vector<TableView<TestTable>::IndexType> evenIndicesTest = {0, 2, 4};

    TableView<TestTable> view(&testTable);
    TableView<TestTable> viewOdd = view.filter([](int i) -> bool { return i % 2; });
    
    std::vector<TableView<TestTable>::IndexType> oddIndices;
    oddIndices.insert(oddIndices.end(), viewOdd.begin(), viewOdd.end());

    compareTwoVectors(oddIndices, oddIndicesTest);

    std::vector<TableView<TestTable>> views = view.split([](int i) -> int { return i % 2; }, 2);
    EXPECT_EQ(views.size(), 2);

    oddIndices.clear();
    oddIndices.insert(oddIndices.end(), views.at(1).begin(), views.at(1).end());

    compareTwoVectors(oddIndices, oddIndicesTest);

    std::vector<TableView<TestTable>::IndexType> evenIndices;
    evenIndices.insert(evenIndices.end(), views.at(0).begin(), views.at(0).end());

    compareTwoVectors(evenIndices, evenIndicesTest);
}