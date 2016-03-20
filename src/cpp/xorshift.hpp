// http://gamedev.stackexchange.com/questions/25730/random-numbers-on-c

#include <limits>

// Simple shift-and-xor random number generator. This
// implementation provides an optimal period for its storage size,
// passes most empirical tests, and is faster and smaller than
// more popular approaches like MT.
//
// However, it is insecure. Don't use it as a source of
// cryptographic randomness.
//
// See www.jstatsoft.org/v08/i14/paper for the algorithm, and
// www.open-std.org/jtc1/sc22/wg21/docs/papers/2003/n1452.html
// for the structure.
class xorshift {

public:
    typedef std::uint32_t result_type;

    struct state_type {
        result_type x;
        result_type y;
        result_type z;
        result_type w;
    };

    xorshift(void)
			: state_({123456789, 362436069, 521288629, 88675123}) {}
    explicit xorshift(result_type r)
			: state_({123456789, 362436069, 521288629, r}) {}
    explicit xorshift(const state_type &seed)
			: state_(seed) {}

    void seed(void) {
			state({123456789, 362436069, 521288629, 88675123});
		}

    void seed(result_type r) {
			state({123456789, 362436069, 521288629, r});
		}

    void seed(const state_type &s) {
			state(s);
		}

		// Generate a uniformly-distributed random integer of
		// result_type.
    result_type operator()(void) {
			result_type t = state_.x ^ (state_.x << 15);
			state_.x = state_.y; state_.y = state_.z; state_.z = state_.w;
			return state_.w = state_.w ^ (state_.w >> 21) ^ (t ^ (t >> 4));
		}

    // Discard the next z random values.
    void discard(unsigned long long z) {
			while (z--) (*this)();
		}

		// Get or set the entire state. This can be used to store
		// and later re-load the state in any format.
    const state_type &state(void) const { return state_; }
    void state(const state_type &state) { state_ = state; }

    // Random number bounds. Used by random distributions; you
    // probably don't need to call these directly.
    static result_type min(void) { return std::numeric_limits<result_type>::min(); }
    static result_type max(void) { return std::numeric_limits<result_type>::max(); }

private:
    state_type state_;
};

//// Two engines compare as equal if their states are
//// bitwise-identical, i.e. if they would generate the same
//// numbers forever.
//bool operator==(const xorshift &lhs, const xorshift &rhs);
//bool operator!=(const xorshift &lhs, const xorshift &rhs);
