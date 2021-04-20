#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <shapelet_transform.hpp>

//TODO: add tests...
TEST_CASE( "Quick check", "[main]" ) {
    REQUIRE( 4.5 == Approx(4.666666) );
}