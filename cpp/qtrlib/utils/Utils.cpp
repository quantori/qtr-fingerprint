#include "Utils.h"

namespace qtr {

bool endsWith(const std::string &a, const std::string &b) {
    return a.size() >= b.size() &&
           (0 == a.compare(a.size() - b.size(), b.size(), b));
}

void emptyArgument(const std::string &argument, const std::string &message) {
    if (argument.empty()) {
        LOG(ERROR) << message;
        exit(-1);
    }
}

int chexToInt(char letter) {
    return letter >= '0' && letter <= '9' ? letter - '0' : letter - 'a' + 10;
}

} // namespace qtr