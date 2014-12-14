#include <cassert>

#include <complex>
#include <algorithm>

#include "sdft/sdft.h"

//
// Exported state struct and internal inheriting struct definitions, templated on the floating point type.
//

struct sdft_State {
    virtual enum sdft_Error validate() = 0;

    virtual enum sdft_Error push_next_sample(void *next_sample) = 0;

    virtual void *get_spectrum() = 0;

    virtual void *unshift_and_get_signal() = 0;

    virtual enum sdft_Error combine_with(struct sdft_State *other, void *buffer) = 0;

    virtual ~sdft_State()
    {
    };
};

template<typename Float>
struct Combined;

template<typename Float>
struct Impl : public sdft_State {
    Impl(void *signal, void *spectrum, void *phase_offsets, size_t number_of_samples,
            enum sdft_SignalTraits signal_traits);

    sdft_Error validate();

    void clear();

    size_t get_number_of_samples() const
    {
        return _number_of_samples;
    }

    sdft_Error push_next_sample(void *next_sample);

    void *get_spectrum()
    {
        return _spectrum;
    }

    void *unshift_and_get_signal();

    virtual sdft_Error combine_with(struct sdft_State *other, void *buffer)
    {
        Impl<Float> *o = dynamic_cast<Impl<Float> *>(other);
        if (o != 0 && o->_number_of_samples == _number_of_samples && o->_signal_traits == _signal_traits) {
            new(buffer) Combined<Float>(this, o);
            return SDFT_NO_ERROR;
        }

        return SDFT_NOT_COMBINABLE;
    }

private:
    typedef std::complex<Float> cplx;

    bool matches_signal_trait(const cplx &c) const;

    cplx *_signal;
    cplx *_spectrum;
    cplx *_phase_offsets;
    size_t _signal_index;
    size_t _number_of_samples;
    enum sdft_SignalTraits _signal_traits;
};

template<typename Float>
struct Combined : public sdft_State {
    Combined(Impl<Float> *first, Impl<Float> *second);

    sdft_Error validate();

    sdft_Error push_next_sample(void *next_sample);

    void *unshift_and_get_signal();

    void *get_spectrum()
    {
        assert(_clear_counter <= 2 * _number_of_samples);
        // See the invariant in push_next_sample.
        return _clear_counter <= _number_of_samples
                ? _first->get_spectrum()
                : _second->get_spectrum();
    }

    virtual enum sdft_Error combine_with(struct sdft_State *other, void *buffer)
    {
        return SDFT_NOT_COMBINABLE;
    }

private:
    Impl<Float> *_first;
    Impl<Float> *_second;
    size_t _number_of_samples;
    size_t _clear_counter;
};

//
// Implementations of exported functions
//

size_t sdft_size_of_state()
{
    return std::max(sizeof(struct Impl<long double>), sizeof(struct Combined<long double>));
}

enum sdft_Error sdft_init_from_buffers(
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
            new(s) Impl<float>(signal, spectrum, buffer, number_of_samples, signal_traits);
            break;
        case SDFT_DOUBLE:
            new(s) Impl<double>(signal, spectrum, buffer, number_of_samples, signal_traits);
            break;
        case SDFT_LONG_DOUBLE:
            new(s) Impl<long double>(signal, spectrum, buffer, number_of_samples, signal_traits);
            break;
    }

    return s->validate();
}

enum sdft_Error sdft_init_combine(struct sdft_State *state, struct sdft_State *first, struct sdft_State *second)
{
    sdft_Error err = first->combine_with(second, state);
    if (err != SDFT_NO_ERROR) {
        return err;
    }

    return state->validate();
}

enum sdft_Error sdft_push_next_sample(struct sdft_State *s, void *next_sample)
{
    return s->push_next_sample(next_sample);
}

void *sdft_get_spectrum(struct sdft_State *s)
{
    return s->get_spectrum();
}

void *sdft_unshift_and_get_signal(struct sdft_State *s)
{
    return s->unshift_and_get_signal();
}

//
// Templated implementations of the precision dependent functions
//

template<typename Float>
Impl<Float>::Impl(void *signal, void *spectrum, void *phase_offsets,
        size_t number_of_samples, enum sdft_SignalTraits signal_traits)
        : _signal((cplx *) signal), _spectrum((cplx *) spectrum), _phase_offsets((cplx *) phase_offsets),
          _signal_index(0), _number_of_samples(number_of_samples), _signal_traits(signal_traits)
{
    // generate the phase offsets
    const Float double_pi = static_cast<Float>(2 * 3.141592653589793238462643383279502884); // Enough precision for everyone!
    for (size_t i = 0; i < number_of_samples; i++) {
        cplx angle(0, double_pi * i / number_of_samples);
        _phase_offsets[i] = exp(angle);
    };
}

template<typename Float>
bool Impl<Float>::matches_signal_trait(typename Impl::cplx const &c) const
{
    return _signal_traits == SDFT_REAL_AND_IMAG
            || (_signal_traits == SDFT_REAL_ONLY && imag(c) == 0)
            || (_signal_traits == SDFT_IMAG_ONLY && real(c) == 0);
}

template<typename Float>
sdft_Error Impl<Float>::validate()
{
    // the window should be of at least length 1
    if (_number_of_samples < 1) {
        return SDFT_WINDOW_TOO_SHORT;
    }

    // check for violations of the signal traits
    for (size_t i = 0; i < _number_of_samples; ++i) {
        if (!matches_signal_trait(_signal[i])) {
            return SDFT_SIGNAL_TRAIT_VIOLATION;
        }
    }

    return SDFT_NO_ERROR;
}

template<typename Float>
void Impl<Float>::clear()
{
    for (size_t i = 0; i < _number_of_samples; ++i) {
        _signal[i] = 0;
    }

    for (size_t i = 0; i < _number_of_samples; ++i) {
        _spectrum[i] = 0;
    }

    _signal_index = 0;
}

template<typename Float>
sdft_Error Impl<Float>::push_next_sample(void *next_sample)
{
    assert(_signal_index >= 0 && _signal_index < _number_of_samples);

    cplx ns = *(cplx *) next_sample;

    if (!matches_signal_trait(ns)) {
        return SDFT_SIGNAL_TRAIT_VIOLATION;
    }

    cplx delta = ns - _signal[_signal_index];

    size_t n_samples = _number_of_samples;
    size_t n_bins = _signal_traits == SDFT_REAL_AND_IMAG
            ? _number_of_samples
            : _number_of_samples / 2; // only first half of spectrum relevant

    for (size_t i = 0; i < n_bins; ++i) {
        _spectrum[i] = (_spectrum[i] + delta) * _phase_offsets[i];
    }

    _signal[_signal_index] += delta;
    if (++_signal_index == n_samples) {
        _signal_index = 0;
    }

    return SDFT_NO_ERROR;
}

template<typename Float>
void *Impl<Float>::unshift_and_get_signal()
{
    size_t n = _number_of_samples;

    if (_signal_index == 0) {
        // we are lucky, there's nothing to shift
        return _signal;
    }

    // instead of taking the signal index as offset, we take n - signal_index
    // as offset to follow a cycle through the array.
    assert(n > _signal_index);
    size_t ofs = _signal_index;

    cplx first = _signal[0];
    size_t cur = 0;
    for (size_t i = 0; i < n - 1; ++i) {
        size_t next = (cur + ofs) % n;
        _signal[cur] = _signal[next];
        cur = next;
    }

    assert((cur + ofs) % n == 0); // or: next == 0
    _signal[cur] = first;
	_signal_index = 0;

    return _signal;
}

template<typename Float>
Combined<Float>::Combined(Impl<Float> *first, Impl<Float> *second)
        : _first(first), _second(second), _number_of_samples(first->get_number_of_samples()), _clear_counter(0)
{
    _second->clear();
}

template<typename Float>
sdft_Error Combined<Float>::validate()
{
    return SDFT_NO_ERROR;
}

template<typename Float>
sdft_Error Combined<Float>::push_next_sample(void *next_sample)
{
    // Invariant: 0 <= _clear_counter <= _number_of_samples
    //                  iff _first has the valid spectrum
    //            _number_of_samples < _clear_counter <= 2*_number_of_samples
    //                  iff _second has the valid spectrum
    assert(_clear_counter <= 2 * _number_of_samples);

    if (_clear_counter == _number_of_samples) {
        _first->clear();
    } else if (_clear_counter == 2 * _number_of_samples) {
        _second->clear();
        _clear_counter = 0;
    }

    enum sdft_Error err = _first->push_next_sample(next_sample);
    if (err != SDFT_NO_ERROR) {
        return err;
    }

    err = _second->push_next_sample(next_sample);
    if (err != SDFT_NO_ERROR) {
        return err;
    }

    _clear_counter++;

    return SDFT_NO_ERROR;
}

template<typename Float>
void *Combined<Float>::unshift_and_get_signal()
{
    assert(_clear_counter <= 2 * _number_of_samples);
    // See the invariant in push_next_sample.
    return _clear_counter <= _number_of_samples
            ? _first->unshift_and_get_signal()
            : _second->unshift_and_get_signal();
}
