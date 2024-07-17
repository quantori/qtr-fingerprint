#pragma once

#include "indigo.h"
#include "base_cpp/array.h"

CEXPORT int qtrBingoSearchSub(int db, int query_obj, const char *options, const indigo::Array <byte> &fp);

CEXPORT int qtrAromatize(int db, int query_obj);