#include "SplitterTree.h"

using namespace std;
using namespace qtr;

uint64_t qtr::SplitterTree::_countOfNodes = 0;


/**
 * Call func for each record in data file
 * @tparam Functor call functor on IndigoFingerprint and string that represent smiles
 * @param filename
 * @param readSmiles if false than smiles equals to an empty string
 * @param func
 */
template<typename Functor>
void forEachLine(const string &filename, bool readSmiles, Functor &&func) {
    ifstream input(filename);
    while (!input.eof()) { //TODO parallelize reading?
        IndigoFingerprint currFP;
        // Read fingerprint
        for (uint64_t i = 0, j = 0;
             i < currFP.sizeInBytes; ++i, j += BIT_IN_BYTE) { // TODO Read fingerprint with one read
            auto curr = input.get();
            for (uint64_t k = 0; k < BIT_IN_BYTE; ++k)
                currFP[j + k] = (curr & (1ull << k));
        }
        // Read smiles
        if (readSmiles) { // TODO buffer reading
            string smiles = "";
            char c;
            while ((c = input.get()) != '\n')
                smiles += c;
            func(currFP, smiles);
        } else {
            while (input.get() != '\n') {}
            func(currFP, "");
        }
    }
}

vector<string> SplitterTree::split(size_t maxDepth, uint64_t minCountInNode) {
    if (_depth > maxDepth)
        return {_filename};

    // Finding the best column to split
    uint64_t countOfRecords = 0;
    size_t countOfBitsInFP = fromBytesToBits(IndigoFingerprint::sizeInBytes);
    uint64_t onesInColumn[countOfBitsInFP];
    memset(onesInColumn, 0, sizeof(onesInColumn));
    forEachLine(_filename, false,
                [&countOfRecords, &countOfBitsInFP, &onesInColumn](const IndigoFingerprint &fp, const string &_) {
                    ++countOfRecords;
                    for (size_t i = 0; i < countOfBitsInFP; ++i)
                        onesInColumn[i] += fp[i];
                });

    if (countOfRecords < minCountInNode) // TODO check it faster?
        return {_filename};

    uint64_t bestBitSplit = -1;
    uint64_t bestBitSplitDeviation = UINT64_MAX;
    auto countInARangeOfColumns = [&onesInColumn, &countOfRecords](size_t left, size_t right)
            -> pair<uint64_t, uint64_t> {
        uint64_t bestBitSplit = -1;
        uint64_t bestBitSplitDeviation = UINT64_MAX;
        for (size_t i = left; i < right; ++i) {
            uint64_t deviation = ((2ll * onesInColumn[i]) >= countOfRecords) ?
                                 ((2ll * onesInColumn[i]) - countOfRecords) :
                                 (countOfRecords - (2ll * onesInColumn[i]));
            if (deviation < bestBitSplitDeviation) {
                bestBitSplitDeviation = deviation;
                bestBitSplit = i;
            }
        }
        return {bestBitSplit, bestBitSplitDeviation};
    };
    if (_depth < COUNT_THREADS_POW) { // Async calculate in each range // TODO extract to func
        size_t availableThreadsCount = (1ull << (COUNT_THREADS_POW - _depth));
        vector<future<pair<uint64_t, uint64_t>>> answers(availableThreadsCount);
        size_t lengthForThread = countOfBitsInFP / availableThreadsCount; // Length of each range
        for (size_t thread = 0; thread < availableThreadsCount; ++thread) {
            size_t left = thread * lengthForThread;
            size_t right = (thread + 1 != availableThreadsCount) ? (left + lengthForThread)
                                                                 : countOfBitsInFP; // if it's last thread, then right should be max
            answers[thread] = async(launch::async, countInARangeOfColumns, left, right);
        }
        for (auto &ans: answers) {
            auto[bit, deviation] = ans.get();
            if (deviation < bestBitSplitDeviation) {
                bestBitSplitDeviation = deviation;
                bestBitSplit = bit;
            }
        }
    } else {
        tie(bestBitSplit, bestBitSplitDeviation) = countInARangeOfColumns(0, countOfBitsInFP);
    }
    _splitBit = bestBitSplit;

    // create two files for left and right child
    _leftChild = new SplitterTree(to_string(_countOfNodes), _depth + 1);
    _rightChild = new SplitterTree(to_string(_countOfNodes), _depth + 1);
    ofstream outLeft(_leftChild->_filename);
    ofstream outRight(_rightChild->_filename);
    forEachLine(_filename, true,
                [&outLeft, &outRight, &bestBitSplit](const IndigoFingerprint &fp, const string &smiles) {
                    ofstream *writeTo = &outLeft;
                    if (fp[bestBitSplit])
                        writeTo = &outRight;
                    auto toWrite = fp.getBytes();
                    char buffer[toWrite.size()];
                    for (size_t i = 0; i < toWrite.size(); ++i) {
                        buffer[i] = (char) toWrite[i];
                    } // TODO test just buffer = (char*)fp.getBytes()
                    writeTo->write(buffer, toWrite.size())
                            .write(smiles.c_str(), smiles.size())
                            .write("\n", 1);
                });
    vector<string> leftSplitted;
    vector<string> rightSplited;
    if (_depth < COUNT_THREADS_POW) { // Run splitting in children async
        auto launchAsync = [&maxDepth, &minCountInNode](SplitterTree *node) {
            return async(launch::async, &SplitterTree::split,
                         node,
                         maxDepth, minCountInNode);
        };
        auto futureLeftSplitted = launchAsync(_leftChild);
        auto futureRightSplitted = launchAsync(_rightChild);
        leftSplitted = std::move(futureLeftSplitted.get());
        rightSplited = std::move(futureRightSplitted.get());
    } else { // Not async
        leftSplitted = std::move(_leftChild->split(maxDepth, minCountInNode));
        rightSplited = std::move(_rightChild->split(maxDepth, minCountInNode));
    }
    remove(_filename.c_str()); // Delete current file
    if (leftSplitted.size() < rightSplited.size()) // Adding smaller one to the bigger one
        swap(leftSplitted, rightSplited);
    leftSplitted.insert(leftSplitted.begin(), rightSplited.begin(),
                        rightSplited.end());
    return leftSplitted;
}

uint64_t SplitterTree::getNum(SplitterTree *node) {
    return (node == nullptr) ? (uint64_t) -1 : node->_nodeNumber;
}

void SplitterTree::saveTo(std::ostream &out) {
    out.write((char *) &_nodeNumber, sizeof(_nodeNumber));
    out.write((char *) &_splitBit, sizeof(_splitBit)); // TODO set it to -1 in a leaf?
    uint64_t idLeft = getNum(_leftChild);
    uint64_t idRight = getNum(_rightChild);
    out.write((char *) &idLeft, sizeof(idLeft));
    out.write((char *) &idRight, sizeof(idRight));
    out.write("\n", 1);

    if (_leftChild != nullptr)
        _leftChild->saveTo(out);
    if (_rightChild != nullptr)
        _rightChild->saveTo(out);
}

std::vector<int> SplitterTree::getNonCorrelatingColumns() {
    const size_t countOfColumns = fromBytesToBits(IndigoFingerprint::sizeInBytes);
    double expectedValue[countOfColumns];
    forEachLine(_filename, false, [&expectedValue](const IndigoFingerprint &fp, const string &_) {
        for (size_t i = 0; i < countOfColumns; ++i)
            expectedValue[i] += fp[i];
    });
    for (size_t i = 0; i < countOfColumns; ++i)
        expectedValue[i] /= countOfColumns;

    double correlation[countOfColumns][countOfColumns];
    for (size_t i = 0; i < countOfColumns; ++i) { // TODO parallelize? maybe if we have small amount of leaves [1]
        for (size_t j = i + 1; j < countOfColumns; ++j) {
            correlation[i][j] =
                    abs(1.0 - 2.0 * expectedValue[j] - 2.0 * expectedValue[i] +
                        4.0 * expectedValue[i] * expectedValue[j]); // TODO fix formula
            correlation[j][i] = correlation[i][j];
        }
        correlation[i][i] = 0;
    }
    vector<pair<double, int>> result(countOfColumns); // Correlation coefficient and it's column number
    for (size_t i = 0; i < countOfColumns; ++i) // TODO see [1]
        result[i] = {*max_element(correlation[i], correlation[i] + countOfColumns), i};
    sort(result.begin(), result.end());
    vector<int> numbers(countOfColumns);
    for (size_t i = 0; i < countOfColumns; ++i)
        numbers[i] = result[i].second;
    return numbers;
}

void SplitterTree::saveNonCorrelatingColumnInEachBucket() {
    if (_leftChild == nullptr && _rightChild == nullptr) { // it's a leaf
        ofstream out(_filename + "OrderColumns");
        auto columnOrder = getNonCorrelatingColumns();
        for (auto &column: columnOrder)
            out << column << ' ';
    } else {
        if (_depth < COUNT_THREADS_POW) { // we have free threads
            future<void> waitLeft;
            future<void> waitRight;
            if (_leftChild != nullptr)
                waitLeft = async(launch::async, &SplitterTree::saveNonCorrelatingColumnInEachBucket, _leftChild);
            if (_rightChild != nullptr)
                waitRight = async(launch::async, &SplitterTree::saveNonCorrelatingColumnInEachBucket, _rightChild);
            waitLeft.get();
            waitRight.get();
        } else {
            if (_leftChild != nullptr)
                _leftChild->saveNonCorrelatingColumnInEachBucket();
            if (_rightChild != nullptr)
                _rightChild->saveNonCorrelatingColumnInEachBucket();
        }
    }
}
