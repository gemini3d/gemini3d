#include <cstdint>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>

#include <boost/ut.hpp>

// based on template from example/matcher.cpp of Boost.UT
// https://github.com/boost-ext/ut/blob/master/example/matcher.cpp
// but modified considerably to duck-type string-like and integer-like types

template <class... Ts>
class any_of {
public:
    constexpr explicit any_of(Ts... ts) : ts_{ts...} {}

    constexpr auto operator==(std::string_view t) const {
        return std::apply([t](const auto&... args) {
            return eq(t, to_sv(args)...);
        }, ts_);
    }

    template <class T,
              class = std::enable_if_t<
                  std::is_integral_v<std::remove_cv_t<std::remove_reference_t<T>>> ||
                  std::is_enum_v<std::remove_cv_t<std::remove_reference_t<T>>>
              >>
    constexpr auto operator==(T t) const {
        return std::apply([t](const auto&... args) {
            return eq_integer(as_integer(t), as_integer(args)...);
        }, ts_);
    }

private:
    template <class T>
    using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

    template <class T>
    static constexpr bool is_integer_like_v =
        std::is_integral_v<remove_cvref_t<T>> || std::is_enum_v<remove_cvref_t<T>>;

    template <class T>
    static constexpr auto as_integer(T value) {
        using value_t = remove_cvref_t<T>;
        if constexpr (std::is_enum_v<value_t>) {
            return static_cast<std::underlying_type_t<value_t>>(value);
        } else {
            return value;
        }
    }

    static constexpr std::string_view to_sv(std::string_view value) {
        return value;
    }

    static std::string_view to_sv(const std::string& value) {
        return value;
    }

    static std::string_view to_sv(const char* value) {
        return value ? std::string_view{value} : std::string_view{};
    }

    template <class T, class U, class... TArgs>
    static constexpr auto eq(const T& t, const U& u, const TArgs&... args) {
        using namespace boost::ut;
        if constexpr (sizeof...(args) > 0) {
            return (that % detail::value{u} == t) or eq(t, args...);
        } else {
            return (that % detail::value{u} == t);
        }
    }

    template <class T, class U>
    static constexpr bool integer_equal(T lhs, U rhs) {
        static_assert(is_integer_like_v<T> && is_integer_like_v<U>, "integer_equal expects integer-like types");

        if constexpr (std::is_signed_v<T> == std::is_signed_v<U>) {
            return lhs == rhs;
        } else if constexpr (std::is_signed_v<T>) {
            if (lhs < 0) {
                return false;
            }
            return static_cast<std::uintmax_t>(lhs) == static_cast<std::uintmax_t>(rhs);
        } else {
            if (rhs < 0) {
                return false;
            }
            return static_cast<std::uintmax_t>(lhs) == static_cast<std::uintmax_t>(rhs);
        }
    }

    template <class T, class U, class... TArgs>
    static constexpr bool eq_integer(const T& t, const U& u, const TArgs&... args) {
        if constexpr (sizeof...(args) > 0) {
            return integer_equal(t, u) || eq_integer(t, args...);
        } else {
            return integer_equal(t, u);
        }
    }

    std::tuple<Ts...> ts_;
};
