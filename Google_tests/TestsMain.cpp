//
// Created by florian on 24.10.23.
//


// AllTests.cpp
#include "gtest/gtest.h"
#include "GraphClassTest.h"
#include "GraphAlgorithmsTests.h"
#include "GraphClosureTests.h"

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}