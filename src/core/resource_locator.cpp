/**
 * @file resource_locator.cpp
 * @brief Implementation of the bundled-resources locator.
 */
#include <pairfinder/core/resource_locator.hpp>

#include <cstdlib>
#include <vector>

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#elif defined(__linux__)
#include <system_error>
#endif

namespace pairfinder::resources {

namespace {

// Absolute path to the running executable (empty if unavailable).
std::filesystem::path executable_path() {
#if defined(__APPLE__)
    std::uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);  // first call reports needed size
    std::string buf(size, '\0');
    if (_NSGetExecutablePath(buf.data(), &size) != 0) return {};
    buf.resize(std::char_traits<char>::length(buf.c_str()));
    std::error_code ec;
    auto canon = std::filesystem::weakly_canonical(std::filesystem::path(buf), ec);
    return ec ? std::filesystem::path(buf) : canon;
#elif defined(__linux__)
    std::error_code ec;
    auto p = std::filesystem::read_symlink("/proc/self/exe", ec);
    return ec ? std::filesystem::path{} : p;
#else
    return {};
#endif
}

// A directory is a valid resources root if it carries our sentinels.
bool valid_root(const std::filesystem::path& p) {
    if (p.empty()) return false;
    std::error_code ec;
    return std::filesystem::exists(p / "templates", ec) &&
           std::filesystem::exists(p / "richardson_suites.json", ec);
}

std::filesystem::path resolve_root() {
    if (const char* env = std::getenv("PARSE_RESOURCES")) {
        std::filesystem::path p(env);
        if (valid_root(p)) return p;
    }
    const std::filesystem::path exe = executable_path();
    if (!exe.empty()) {
        const std::filesystem::path dir = exe.parent_path();
        for (const char* rel : {"resources", "../resources", "../share/parse/resources"}) {
            std::filesystem::path p = dir / rel;
            if (valid_root(p)) {
                std::error_code ec;
                auto canon = std::filesystem::weakly_canonical(p, ec);
                return ec ? p : canon;
            }
        }
    }
    for (const char* rel : {"resources", "../resources", "../../resources"}) {
        std::filesystem::path p(rel);
        if (valid_root(p)) return std::filesystem::absolute(p);
    }
#ifdef PAIRFINDER_SOURCE_RESOURCES
    {
        std::filesystem::path p(PAIRFINDER_SOURCE_RESOURCES);
        if (valid_root(p)) return p;
    }
#endif
    return {};
}

}  // namespace

const std::filesystem::path& root() {
    static const std::filesystem::path r = resolve_root();
    return r;
}

std::filesystem::path file(const std::string& relative) { return root() / relative; }
std::filesystem::path templates_dir() { return root() / "templates"; }
std::filesystem::path config(const std::string& name) { return root() / "config" / name; }
std::filesystem::path reference_dir() { return root() / "reference"; }
std::filesystem::path basepair_idealized() { return root() / "basepair-idealized"; }
std::filesystem::path basepair_exemplars() { return root() / "basepair-catalog-exemplars"; }

}  // namespace pairfinder::resources
