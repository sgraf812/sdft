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

    virtual void *unshift_and_get_window() = 0;

    virtual enum sdft_Error combine_with(struct sdft_State *other, void *buffer) = 0;

    virtual ~sdft_State()
    {
    };
};

template<typename Float>
struct Combined;

template<typename Float>
struct Impl : public sdft_State {
    Impl(void *signal, void *spectrum, void *phase_offsets, size_t window_size,
            enum sdft_SignalTraits signal_traits);

    sdft_Error validate();

    void clear();

    size_t get_window_size() const
    {
        return _window_size;
    }

    sdft_Error push_next_sample(void *next_sample);

    void *get_spectrum()
    {
        return _spectrum;
    }

    void *unshift_and_get_window();

    sdft_Error combine_with(struct sdft_State *other, void *buffer)
    {
        Impl<Float> *o = dynamic_cast<Impl<Float> *>(other);
        if (o != 0 && o->_window_size == _window_size && o->_signal_traits == _signal_traits) {
            new(buffer) Combined<Float>(this, o);
            return SDFT_NO_ERROR;
        }

        return SDFT_NOT_COMBINABLE;
    }

private:
    typedef std::complex<Float> cplx;

    bool matches_signal_trait(const cplx &c) const;

    cplx *_window;
    cplx *_spectrum;
    cplx *_phase_offsets;
    size_t _window_index;
    size_t _window_size;
    enum sdft_SignalTraits _signal_traits;
};

template<typename Float>
struct Combined : public sdft_State {
    Combined(Impl<Float> *first, Impl<Float> *second);

    sdft_Error validate();

    sdft_Error push_next_sample(void *next_sample);

    void *unshift_and_get_window();

    void *get_spectrum()
    {
        assert(_clear_counter <= 2 * _window_size);
        // See the invariant in push_next_sample.
        return _clear_counter <= _window_size
                ? _first->get_spectrum()
                : _second->get_spectrum();
    }

    enum sdft_Error combine_with(struct sdft_State *, void *)
    {
        return SDFT_NOT_COMBINABLE;
    }

private:
    Impl<Float> *_first;
    Impl<Float> *_second;
    size_t _window_size;
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
        void *window,
        void *spectrum,
        void *buffer,
        size_t window_size,
        enum sdft_SignalTraits signal_traits)
{

    switch (precision) {
        case SDFT_SINGLE:
            new(s) Impl<float>(window, spectrum, buffer, window_size, signal_traits);
            break;
        case SDFT_DOUBLE:
            new(s) Impl<double>(window, spectrum, buffer, window_size, signal_traits);
            break;
        case SDFT_LONG_DOUBLE:
            new(s) Impl<long double>(window, spectrum, buffer, window_size, signal_traits);
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

void *sdft_unshift_and_get_window(struct sdft_State *s)
{
    return s->unshift_and_get_window();
}

//
// Templated implementations of the precision dependent functions
//

template<typename Float>
Impl<Float>::Impl(void *signal, void *spectrum, void *phase_offsets,
        size_t window_size, enum sdft_SignalTraits signal_traits)
        : _window((cplx *) signal), _spectrum((cplx *) spectrum), _phase_offsets((cplx *) phase_offsets),
          _window_index(0), _window_size(window_size), _signal_traits(signal_traits)
{
    // generate the phase offsets
    const Float double_pi = static_cast<Float>(2 * 3.141592653589793238462643383279502884); // Enough precision for everyone!
    for (size_t i = 0; i < window_size; i++) {
        cplx angle(0, double_pi * i / window_size);
        _phase_offsets[i] = std::exp(angle);
    };
}

template<typename Float>
bool Impl<Float>::matches_signal_trait(typename Impl::cplx const &c) const
{
    return _signal_traits == SDFT_REAL_AND_IMAG
            || (_signal_traits == SDFT_REAL_ONLY && std::imag(c) == 0)
            || (_signal_traits == SDFT_IMAG_ONLY && std::real(c) == 0);
}

template<typename Float>
sdft_Error Impl<Float>::validate()
{
    // the window should be of at least length 1
    if (_window_size < 1) {
        return SDFT_WINDOW_TOO_SHORT;
    }

    // check for violations of the signal traits
    for (size_t i = 0; i < _window_size; ++i) {
        if (!matches_signal_trait(_window[i])) {
            return SDFT_SIGNAL_TRAIT_VIOLATION;
        }
    }

    return SDFT_NO_ERROR;
}

template<typename Float>
void Impl<Float>::clear()
{
    for (size_t i = 0; i < _window_size; ++i) {
        _window[i] = 0;
    }

    for (size_t i = 0; i < _window_size; ++i) {
        _spectrum[i] = 0;
    }

    _window_index = 0;
}

template<typename Float>
sdft_Error Impl<Float>::push_next_sample(void *next_sample)
{
    assert(_window_index >= 0 && _window_index < _window_size);

    cplx ns = *(cplx *) next_sample;

    if (!matches_signal_trait(ns)) {
        return SDFT_SIGNAL_TRAIT_VIOLATION;
    }

    cplx delta = ns - _window[_window_index];

    size_t n_bins = _signal_traits == SDFT_REAL_AND_IMAG
            ? _window_size
            : _window_size / 2; // only first half of spectrum relevant

    for (size_t i = 0; i < n_bins; ++i) {
        _spectrum[i] = (_spectrum[i] + delta) * _phase_offsets[i];
    }

    _window[_window_index] = ns;
    if (++_window_index == _window_size) {
        _window_index = 0;
    }

    return SDFT_NO_ERROR;
}

template<typename Float>
void *Impl<Float>::unshift_and_get_window()
{
    assert(_window_size > _window_index);

    // this cyclically shifts the element at _window_index to the front
    std::rotate(_window, _window + _window_index, _window + _window_size);

    _window_index = 0;
    return _window;
}

template<typename Float>
Combined<Float>::Combined(Impl<Float> *first, Impl<Float> *second)
        : _first(first), _second(second), _window_size(first->get_window_size()), _clear_counter(0)
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
    // Invariant: 0 <= _clear_counter <= _window_size
    //                  iff _first has the valid spectrum
    //            _window_size < _clear_counter <= 2*_window_size
    //                  iff _second has the valid spectrum
    assert(_clear_counter <= 2 * _window_size);

    if (_clear_counter == _window_size) {
        _first->clear();
    } else if (_clear_counter == 2 * _window_size) {
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
void *Combined<Float>::unshift_and_get_window()
{
    assert(_clear_counter <= 2 * _window_size);
    // See the invariant in push_next_sample.
    return _clear_counter <= _window_size
            ? _first->unshift_and_get_window()
            : _second->unshift_and_get_window();
}
