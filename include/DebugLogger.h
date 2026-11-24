#ifndef DEBUG_LOGGER_H
#define DEBUG_LOGGER_H

#include <iostream>
#include <string>

#ifdef DEBUG_MODE
#define DEBUG_LOG(msg) std::cout << "[DEBUG] " << msg << std::endl
#define DEBUG_CHECK_TEMPERATURE(T, context) \
    if (T < 200.0 || T > 4000.0) { \
        std::cerr << "[WARNING] Extreme temperature " << T << "K in " << context << std::endl; \
    }
#define DEBUG_CHECK_PRESSURE(P, context) \
    if (P < 1000.0 || P > 1e8) { \
        std::cerr << "[WARNING] Extreme pressure " << P/1e5 << "bar in " << context << std::endl; \
    }
#else
#define DEBUG_LOG(msg)
#define DEBUG_CHECK_TEMPERATURE(T, context)
#define DEBUG_CHECK_PRESSURE(P, context)
#endif

#endif
