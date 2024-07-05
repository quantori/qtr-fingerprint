#pragma once

#include "BingoNoSQL.h"
#include "IndigoSession.h"
#include "indigo.h"

inline const std::shared_ptr<indigo_cpp::IndigoSession> globalIndigoSession = indigo_cpp::IndigoSession::create();