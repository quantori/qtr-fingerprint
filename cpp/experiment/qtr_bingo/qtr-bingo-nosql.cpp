#include "qtr-bingo-nosql.h"

#include <cstdio>
#include <string>

#include "bingo_index.h"
#include "bingo_internal.h"
#include "indigo_internal.h"
#include "indigo_molecule.h"
#include "indigo_reaction.h"

using namespace indigo;
using namespace bingo;

namespace {
    template<class T>
    class BingoPool {
    public:
        bool has(long long id) const {
            return map.count(id) > 0;
        }

        sf::safe_shared_hide_obj<std::unique_ptr<T>> &at(long long id) {
            return map.at(id);
        }

        const sf::safe_shared_hide_obj<std::unique_ptr<T>> &at(long long id) const {
            return map.at(id);
        }

        void insert(long long id, std::unique_ptr<T> &&obj) {
            map[id] = std::move(sf::safe_shared_hide_obj<std::unique_ptr<T>>(std::move(obj)));
        }

        long long insert(std::unique_ptr<T> &&obj) {
            map[next_id] = std::move(sf::safe_shared_hide_obj<std::unique_ptr<T>>(std::move(obj)));
            return next_id++;
        }

        void remove(long long id) {
            map.erase(id);
        }

        long long getNextId() {
            return next_id++;
        }

    private:
        std::unordered_map<long long, sf::safe_shared_hide_obj<std::unique_ptr<T>>> map;
        long long next_id = 1;
    };

    struct SearchesData {
        BingoPool<Matcher> searches;
        std::unordered_map<long long, long long> db;
    };

    static sf::safe_shared_hide_obj<BingoPool<BaseIndex>> &_indexes() {
        static sf::safe_shared_hide_obj<BingoPool<BaseIndex>> indexes;
        return indexes;
    }

    static sf::safe_shared_hide_obj<SearchesData> &_searches_data() {
        static sf::safe_shared_hide_obj<SearchesData> searches_data;
        return searches_data;
    }

}

CEXPORT int qtrAromatize(int db, int query_obj) {
    BINGO_BEGIN_DB(db)
                {
                    auto obj_ptr = std::unique_ptr<IndigoObject>(self.getObject(query_obj).clone());
                    IndigoObject &obj = *obj_ptr;
                    assert(IndigoQueryMolecule::is(obj));
                    obj.getBaseMolecule().aromatize(self.arom_options);
                }
    BINGO_END(-1);
}


CEXPORT int qtrBingoSearchSub(int db, int query_obj, const char *options, const indigo::Array<byte> &fp) {
    BINGO_BEGIN_DB(db)
                {
                    auto obj_ptr = std::unique_ptr<IndigoObject>(self.getObject(query_obj).clone());
                    IndigoObject &obj = *obj_ptr;

                    assert(IndigoQueryMolecule::is(obj));

//                    obj.getBaseMolecule().aromatize(self.arom_options);
                    std::unique_ptr<MoleculeSubstructureQueryData> query_data = std::make_unique<MoleculeSubstructureQueryData>(
                            obj.getQueryMolecule());

                    auto matcher = [&]() {
                        const auto bingo_indexes = sf::slock_safe_ptr(_indexes());
                        const auto bingo_index_ptr = sf::slock_safe_ptr(bingo_indexes->at(db));
                        auto &base_index = **bingo_index_ptr;
                        auto molecule_index = dynamic_cast<MoleculeIndex *>(&base_index);
                        assert(molecule_index != nullptr);
//                        return createMatcher(*molecule_index, "sub", query_data.release(), options, fp); // TODO: pass fingerprint in the matcher
                        return (*bingo_index_ptr)->createMatcher("sub", query_data.release(), options);
                    }();

                    {
                        auto searches_data = sf::xlock_safe_ptr(_searches_data());
                        auto search_id = searches_data->searches.insert(std::move(matcher));
                        searches_data->db[search_id] = db;
                        return search_id;
                    }
                }
    BINGO_END(-1);
}
