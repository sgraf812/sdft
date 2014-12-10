#include <complex.h>

#include "sdft/sdft.h"

typedef float complex cplxf;
typedef double complex cplx;
typedef long double complex cplxl;

struct sdft_State {
    void *signal;
    void *spectrum;
    void *phase_offsets;
    size_t signal_index;
    size_t number_of_samples;
    void (*push_next_sample)(struct sdft_State*, void *next_sample);
    void (*unshift_signal)(struct sdft_State*);
    enum sdft_SignalTraits signal_traits;
};

size_t sdft_size_of_state()
{
    return sizeof(struct sdft_State);
}

static size_t number_of_bins(struct sdft_State *s)
{
    return s->signal_traits == SDFT_REAL_AND_IMAG
            ? s->number_of_samples
            : s->number_of_samples / 2;
}

// We fake C++ templates with macro definitions and inclusions of templated_function.h.
// That way, we can support all three floating point precisions at the cost of a virtual dispatch
// on every function call. The runtime cost of the algorithm should far outhweigh that, however.
#define APPEND_TYPE_SUFFIX(x) x##f
#include "templated_functions.h"
#undef APPEND_TYPE_SUFFIX
#define APPEND_TYPE_SUFFIX(x) x
#include "templated_functions.h"
#undef APPEND_TYPE_SUFFIX
#define APPEND_TYPE_SUFFIX(x) x##l
#include "templated_functions.h"
#undef APPEND_TYPE_SUFFIX

void sdft_init(
        struct sdft_State *s,
        enum sdft_FloatPrecision precision,
        void *signal,
        void *spectrum,
        void *buffer,
        size_t number_of_samples,
        enum sdft_SignalTraits signal_traits)
{
    switch (precision) {
        case SDFT_SINGLE:
            impl_initf(s, signal, spectrum, buffer, number_of_samples, signal_traits);
            break;
        case SDFT_DOUBLE:
            impl_init(s, signal, spectrum, buffer, number_of_samples, signal_traits);
            break;
        case SDFT_LONG_DOUBLE:
            impl_initl(s, signal, spectrum, buffer, number_of_samples, signal_traits);
            break;
    }
}

void sdft_push_next_sample(struct sdft_State *s, void *next_sample)
{
    s->push_next_sample(s, next_sample);
}

void sdft_unshift_signal(struct sdft_State *s)
{
    s->unshift_signal(s);
}
