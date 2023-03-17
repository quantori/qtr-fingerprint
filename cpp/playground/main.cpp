#include "bits/stdc++.h"

#include "IndigoSession.h"
#include "IndigoQueryMolecule.h"
#include "IndigoMolecule.h"
#include "IndigoWriteBuffer.h"
#include "IndigoSDFileIterator.h"
#include "indigo_savers.h"
#include "indigo_internal.h"
#include "gzip/gzip_output.h"
#include "gzip/gzip_scanner.h"

#include "string_table_io/StringTableReader.h"
#include "Utils.h"

using namespace std;
using namespace indigo_cpp;
using namespace qtr;

int main() {

    auto indigoSessionPtr = IndigoSession::create();
    TimeMeasurer timeMeasurer;

    filesystem::path filePath("/home/Vsevolod.Vaskin/qtr-fingerprint/cpp/playground/molecules.sdf");
    filesystem::path stPath("/home/Vsevolod.Vaskin/qtr-fingerprint/data/parsed_data/smilesTables/lib1.st");
    int compressRate = 3;


    vector<pair<uint64_t, IndigoMolecule>> molecules;
    {
        TimeMeasurer::FunctionTimeMeasurer timer(timeMeasurer, "st file reading");

        for (const auto& [id, smiles] : qtr::StringTableReader(stPath)) {
            auto mol = indigoSessionPtr->loadMolecule(smiles);
            mol.aromatize();
            molecules.emplace_back(id, mol);
        }
    }

    {
        TimeMeasurer::FunctionTimeMeasurer timer(timeMeasurer, "gzip file writing");

        FileOutput out(filePath.c_str());
        GZipOutput zipOut(out, compressRate);
        for (auto& [id, mol] : molecules) {
            if (zipOut.tell() != 0)
                zipOut.writeChar('\n');
            mol.aromatize();
//            cerr << id << ' ' << mol.canonicalSmiles() << endl;
            zipOut.printf("%d %s$$$$", id, mol.molfile().c_str());
        }
    }

    vector<pair<uint64_t, IndigoMolecule>> molecules2;
    {
        TimeMeasurer::FunctionTimeMeasurer timer(timeMeasurer, "gzip file reading");

        for (const auto &mol: indigoSessionPtr->iterateSDFile(filePath)) {
            auto id = stoi(mol->name());
            mol->aromatize();
            molecules2.emplace_back(id, *mol);
        }
    }

    {
        TimeMeasurer::FunctionTimeMeasurer timer(timeMeasurer, "checking");

        map<uint64_t, IndigoMolecule> molecules2map(molecules2.begin(), molecules2.end());

        if(molecules.size() != molecules2map.size())
            throw exception();
        for (auto& [id, mol] : molecules) {
            auto it = molecules2map.find(id);
            if (it == molecules2map.end())
                throw exception();

            if (!indigoExactMatch(mol.id(), it->second.id(), "ELE"))
                throw exception();

//            cerr << "ok" << endl;
        }
    }


    cout << "Compress rate " << compressRate << endl;

    for (auto& [label, time] : timeMeasurer) {
        cout << label << ' ' << time << endl;
    }

    return 0;
}