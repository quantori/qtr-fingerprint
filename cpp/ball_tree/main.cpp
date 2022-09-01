#include "ColumnsIO.h"

int main(int argc, char *argv[]) {
    {
        qtr::ColumnsWriter writer("/home/Vsevolod.Vaskin/qtr-fingerprint/cpp/test_cols.txt");
        writer.write(10);
        writer.write(20);
        writer.write(0);
        writer.write({100, 200, 300, 400, 500});
    }
    {
        qtr::ColumnsReader reader("/home/Vsevolod.Vaskin/qtr-fingerprint/cpp/test_cols.txt");
        assert(reader.readOne() == 10);
        assert(reader.readOne() == 20);
        assert(reader.readOne() == 0);
        auto r = reader.readAll();
        std::vector<size_t> rr = {100, 200, 300, 400, 500};
        assert(r == rr);
        std::cout << r[0] << ' ' << r[1] << ' ' << r[2];
    }
}