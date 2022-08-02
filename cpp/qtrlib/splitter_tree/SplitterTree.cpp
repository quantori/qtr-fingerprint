#include "SplitterTree.h"

using namespace std;
using namespace qtr;

std::atomic_uint64_t qtr::SplitterTree::_countOfNodes = 0;

/**
 * Call func for each record in data file
 * @tparam Functor call functor on IndigoFingerprint and string that represent smiles
 * @param filename
 * @param readSmiles if false than smiles equals to an empty string
 * @param func
 */
template<typename Functor>
void forEachLine(const filesystem::path &filePath, bool readSmiles, Functor &&func) {
    ifstream input(filePath);
    if (input.bad()) {
        std::cerr << "smth with " << filePath.string() << " file\n";
        return;
    }
    uint64_t cntRecords;
    input.read((char *) (&cntRecords), sizeof(cntRecords));
    for (uint64_t record = 0; record < cntRecords; ++record) { //TODO parallelize reading?
        IndigoFingerprint currFP;
        // Read fingerprint
        currFP.readFrom(input);
        // Read smiles
        if (readSmiles) { // TODO buffer reading?
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

uint64_t SplitterTree::findBestBitToSplit() {
    uint64_t countOfRecords = 0;
    size_t countOfBitsInFP = fromBytesToBits(IndigoFingerprint::sizeInBytes);
    uint64_t onesInColumn[countOfBitsInFP];
    memset(onesInColumn, 0, sizeof(onesInColumn));
    forEachLine(_dir / _filename, false,
                [&countOfRecords, &countOfBitsInFP, &onesInColumn](const IndigoFingerprint &fp, const string &_) {
                    ++countOfRecords;
                    for (size_t i = 0; i < countOfBitsInFP; ++i)
                        onesInColumn[i] += fp[i];
                });

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
    uint64_t bestBitSplit = -1;
    uint64_t bestBitSplitDeviation = UINT64_MAX;
    if (_depth < COUNT_THREADS_POW) { // Async calculate in each range
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
    return bestBitSplit;
}

std::pair<uint64_t, uint64_t> SplitterTree::prepareFilesForChildren(uint64_t bitSplit) {
    uint64_t sizeLeft = 0;
    uint64_t sizeRight = 0;
    ofstream outLeft(_dir / _leftChild->_filename);
    ofstream outRight(_dir / _rightChild->_filename);
    outLeft.write((char *) (&sizeLeft), sizeof(sizeLeft));
    outRight.write((char *) (&sizeRight), sizeof(sizeRight));
    forEachLine(_dir / _filename, true,
                [&outLeft, &outRight, &bitSplit, &sizeLeft, &sizeRight](const IndigoFingerprint &fp,
                                                                        const string &smiles) {
                    ofstream *writeTo = &outLeft;
                    if (fp[bitSplit]) {
                        writeTo = &outRight;
                        sizeRight++;
                    } else
                        sizeLeft++;
                    fp.saveBytes(*writeTo);
                    (*writeTo) << smiles << '\n';
                });
    outLeft.seekp(0, ios::beg);
    outLeft.write((char *) (&sizeLeft), sizeof(sizeLeft));
    outRight.seekp(0, ios::beg);
    outRight.write((char *) (&sizeRight), sizeof(sizeRight));

    remove((_dir / _filename).c_str()); // Delete current file
    return {sizeLeft, sizeRight};
}

std::vector<std::string>
SplitterTree::splitInChildren(uint64_t sizeLeft, uint64_t sizeRight, size_t maxDepth, uint64_t minCountInNode) {
    vector<string> leftSplitted;
    vector<string> rightSplited;
    auto launch = [maxDepth, minCountInNode](SplitterTree *node, uint64_t cntRecords) {
        if (cntRecords >= minCountInNode)
            return node->split(maxDepth, minCountInNode);
        return vector<string>{node->_filename};
    };
    if (sizeLeft != 0 && sizeRight != 0) {
        if (_depth < COUNT_THREADS_POW) { // Run splitting in children async
            auto futureLeftSplit = async(launch::async, launch, _leftChild, sizeLeft);
            auto futureRightSplit = async(launch::async, launch, _rightChild, sizeRight);
            leftSplitted = std::move(futureLeftSplit.get());
            rightSplited = std::move(futureRightSplit.get());
        } else { // Not async
            leftSplitted = std::move(launch(_leftChild, sizeLeft));
            rightSplited = std::move(launch(_rightChild, sizeRight));
        }
    }
    if (leftSplitted.size() < rightSplited.size()) // Adding smaller one to the bigger one
        swap(leftSplitted, rightSplited);
    leftSplitted.insert(leftSplitted.begin(), rightSplited.begin(),
                        rightSplited.end());
    return leftSplitted;
}

vector<string> SplitterTree::split(size_t maxDepth, uint64_t minCountInNode) {
    if (_depth >= maxDepth)
        return {_dir / _filename};

    uint64_t bestBitSplit = findBestBitToSplit();
    _splitBit = bestBitSplit;

    _leftChild = new SplitterTree(_dir, to_string(_countOfNodes), _depth + 1);
    _rightChild = new SplitterTree(_dir, to_string(_countOfNodes), _depth + 1);
    auto[sizeLeft, sizeRight] = prepareFilesForChildren(bestBitSplit);

    return splitInChildren(sizeLeft, sizeRight, maxDepth, minCountInNode);
}

uint64_t SplitterTree::getNum(SplitterTree *node) {
    return (node == nullptr) ? (uint64_t) -1 : node->_nodeNumber;
}

void SplitterTree::saveTo(std::ostream &out, bool writeSize) {
    if (writeSize) {
        uint64_t sz = size();
        out.write((char *) (&sz), sizeof(sz));
    }
    out.write((char *) &_nodeNumber, sizeof(_nodeNumber));
    out.write((char *) &_splitBit, sizeof(_splitBit));
    uint64_t idLeft = getNum(_leftChild);
    uint64_t idRight = getNum(_rightChild);
    out.write((char *) &idLeft, sizeof(idLeft));
    out.write((char *) &idRight, sizeof(idRight));

    if (_leftChild != nullptr)
        _leftChild->saveTo(out, false);
    if (_rightChild != nullptr)
        _rightChild->saveTo(out, false);
}

uint64_t SplitterTree::size() const {
    uint64_t sz = 1;
    if (_leftChild != nullptr)
        sz += _leftChild->size();
    if (_rightChild != nullptr)
        sz += _rightChild->size();
    return sz;
}