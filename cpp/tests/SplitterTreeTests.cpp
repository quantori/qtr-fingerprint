#include "SplitterTree.h"

#include "Common.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <queue>
#include <vector>
#include <utility>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>

using namespace qtr;

TEST(SplitterTree, Spliting) {
    SplitterTree tree("asd");
    tree.split(1, 10000);
}