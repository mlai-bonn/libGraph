//
// Created by florian on 24.10.23.
//


// AllTests.cpp
#include "gtest/gtest.h"
#include "tests/GraphClassTest.h"
#include "tests/GraphAlgorithmsTests.h"
#include "tests/GraphClosureTests.h"
#include "tests/GraphCoreTests.h"

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}