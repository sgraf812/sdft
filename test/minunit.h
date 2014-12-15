#pragma once

#define MU_ASSERT(message, test) do { if (!(test)) return message; } while (0)
#define MU_RUN_TESTS(test) do { char *message = test(); \
                                if (message) return message; } while (0)
extern int tests_run;
