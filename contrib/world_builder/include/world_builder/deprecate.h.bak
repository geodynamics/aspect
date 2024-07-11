// if deprecated is supported, which is definitely C++14, but some compilers support it before C++14.
// Once the world builder switches to C++14 this can be simplified.
#if defined(__has_cpp_attribute)
#if __has_cpp_attribute(deprecated)
#define DEPRECATED(msg, func) [[deprecated(msg)]] func
#endif
#else
#if defined(__GNUC__) || defined(__clang__)
#define DEPRECATED(msg, func) func __attribute__ ((deprecated(msg)))
#elif defined(_MSC_VER)
#define DEPRECATED(msg, func) __declspec(deprecated(msg)) func
#endif
#endif