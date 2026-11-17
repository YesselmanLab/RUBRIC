/**
 * @file test_helper.hpp
 * @brief Minimal test harness for the RUBRIC tests (no external dependency).
 *
 * Each test file defines its tests with TEST(name) { ... } and ends with
 * RUBRIC_TEST_MAIN(). The program runs every test, prints one line per test and
 * "OK" when all pass, and exits non-zero on the first failing check of each test
 * (the other tests still run). Checks print the values involved, so a failure
 * explains itself without a debugger.
 *
 *   TEST(superpose_recovers_rotation) {
 *       CHECK(result.valid);
 *       CHECK_NEAR(result.rms, 0.0, 1e-9);
 *       CHECK_EQ(name, std::string("cWW"));
 *   }
 *   RUBRIC_TEST_MAIN()
 */
#ifndef RUBRIC_TESTS_TEST_HELPER_HPP
#define RUBRIC_TESTS_TEST_HELPER_HPP

#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rubric_test {

/// Raised by a failing check; caught by the runner, which reports the test as failed.
struct CheckFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

inline std::vector<std::pair<std::string, std::function<void()>>>& registry() {
    static std::vector<std::pair<std::string, std::function<void()>>> tests;
    return tests;
}

struct Registrar {
    Registrar(const char* name, std::function<void()> fn) {
        registry().emplace_back(name, std::move(fn));
    }
};

template <typename T>
std::string show(const T& value) {
    std::ostringstream os;
    os << value;
    return os.str();
}
inline std::string show(const std::string& value) { return "\"" + value + "\""; }
inline std::string show(const char* value) { return show(std::string(value)); }
inline std::string show(bool value) { return value ? "true" : "false"; }

[[noreturn]] inline void fail(const char* file, int line, const std::string& message) {
    std::ostringstream os;
    os << file << ":" << line << ": " << message;
    throw CheckFailure(os.str());
}

/// Run every registered test. Returns the process exit code.
inline int run_all() {
    int failed = 0;
    for (const auto& [name, fn] : registry()) {
        try {
            fn();
            std::cout << "  pass  " << name << "\n";
        } catch (const CheckFailure& e) {
            ++failed;
            std::cout << "  FAIL  " << name << "\n        " << e.what() << "\n";
        } catch (const std::exception& e) {
            ++failed;
            std::cout << "  FAIL  " << name << "\n        unexpected exception: " << e.what() << "\n";
        }
    }
    if (failed) {
        std::cout << failed << " of " << registry().size() << " tests failed\n";
        return 1;
    }
    std::cout << registry().size() << " tests\nOK\n";
    return 0;
}

/// Directory with the test data (tests/), set by CMake.
inline std::filesystem::path tests_dir() {
#ifdef RUBRIC_TESTS_DIR
    return RUBRIC_TESTS_DIR;
#else
    return std::filesystem::path(__FILE__).parent_path();
#endif
}

}  // namespace rubric_test

#define RUBRIC_TEST_CAT2(a, b) a##b
#define RUBRIC_TEST_CAT(a, b) RUBRIC_TEST_CAT2(a, b)

/// Define a test. The body runs as a function; a failing check ends it.
#define TEST(name)                                                              \
    static void RUBRIC_TEST_CAT(test_, name)();                                 \
    static const rubric_test::Registrar RUBRIC_TEST_CAT(registrar_, name)(      \
        #name, &RUBRIC_TEST_CAT(test_, name));                                  \
    static void RUBRIC_TEST_CAT(test_, name)()

/// Check a condition; the message shows the expression.
#define CHECK(cond)                                                             \
    do {                                                                        \
        if (!(cond)) rubric_test::fail(__FILE__, __LINE__, "CHECK(" #cond ") failed"); \
    } while (0)

/// Check a condition with an explanatory message built from values.
#define CHECK_MSG(cond, msg)                                                    \
    do {                                                                        \
        if (!(cond)) {                                                          \
            std::ostringstream rubric_os_;                                      \
            rubric_os_ << "CHECK(" #cond ") failed: " << msg;                   \
            rubric_test::fail(__FILE__, __LINE__, rubric_os_.str());            \
        }                                                                       \
    } while (0)

/// Check equality; the message shows both values.
#define CHECK_EQ(a, b)                                                          \
    do {                                                                        \
        const auto& rubric_a_ = (a);                                            \
        const auto& rubric_b_ = (b);                                            \
        if (!(rubric_a_ == rubric_b_))                                          \
            rubric_test::fail(__FILE__, __LINE__,                               \
                "CHECK_EQ(" #a ", " #b ") failed: " + rubric_test::show(rubric_a_) + \
                " != " + rubric_test::show(rubric_b_));                         \
    } while (0)

/// Check |a - b| <= tol; the message shows both values and the tolerance.
#define CHECK_NEAR(a, b, tol)                                                   \
    do {                                                                        \
        const double rubric_a_ = static_cast<double>(a);                        \
        const double rubric_b_ = static_cast<double>(b);                        \
        if (!(std::fabs(rubric_a_ - rubric_b_) <= (tol)))                       \
            rubric_test::fail(__FILE__, __LINE__,                               \
                "CHECK_NEAR(" #a ", " #b ") failed: " + rubric_test::show(rubric_a_) + \
                " vs " + rubric_test::show(rubric_b_) + " (tol " +              \
                rubric_test::show(static_cast<double>(tol)) + ")");             \
    } while (0)

/// Check that an expression throws.
#define CHECK_THROWS(expr)                                                      \
    do {                                                                        \
        bool rubric_threw_ = false;                                             \
        try { (void)(expr); } catch (...) { rubric_threw_ = true; }            \
        if (!rubric_threw_)                                                     \
            rubric_test::fail(__FILE__, __LINE__, "CHECK_THROWS(" #expr ") did not throw"); \
    } while (0)

#define RUBRIC_TEST_MAIN() \
    int main() { return rubric_test::run_all(); }

#endif  // RUBRIC_TESTS_TEST_HELPER_HPP
