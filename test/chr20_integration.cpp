//
// Created by Marco on 24/05/2021.
//

#define CATCH_CONFIG_MAIN
#include <catch/catch.hpp>

#include <string>
#include <iostream>
#include <filesystem>

#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>

#include "../common.hpp"

using namespace std;

// https://stackoverflow.com/questions/6163611/compare-two-files/37575457
bool compareFiles(const string& p1, const string& p2) {
    ifstream f1(p1, ifstream::binary|ifstream::ate);
    ifstream f2(p2, ifstream::binary|ifstream::ate);

    if (f1.fail() || f2.fail()) {
        return false; //file problem
    }

    if (f1.tellg() != f2.tellg()) {
        return false; //size mismatch
    }

    //seek back to beginning and use equal to compare contents
    f1.seekg(0, ifstream::beg);
    f2.seekg(0, ifstream::beg);
    return equal(istreambuf_iterator<char>(f1.rdbuf()),
                      istreambuf_iterator<char>(),
                      istreambuf_iterator<char>(f2.rdbuf()));
}


TEST_CASE("chr20 correctness")
{
    SECTION("index")
    {
        REQUIRE(filesystem::exists("../../test/kmc.out.kmc_pre"));
        REQUIRE(filesystem::exists("../../test/kmc.out.kmc_suf"));
        char *argv[] = {"index", "-k", "35", "-r", "43", "-b", "1", "-f", "EUR_AF",
                        "../../test/chr20.fa", "../../test/chr20.vcf", "../../test/kmc.out",
                        NULL};
        int argc = 12;
        int result = index_main(argc, argv);
        REQUIRE(result == 0);
    }

    SECTION("call")
    {
        char *argv[] = {"call", "-k", "35", "-r", "43", "-b", "1", "-f", "EUR_AF",
                        "../../test/chr20.fa", "../../test/chr20.vcf", "../../test/kmc.out",
                        NULL};
        int argc = 12;
        // Vicious hack...
        ofstream outf("out.vcf");
        streambuf *coutbuf = cout.rdbuf();
        cout.rdbuf(outf.rdbuf());
        int result = call_main(argc, argv);
        cout.flush();
        cout.rdbuf(coutbuf);
        REQUIRE(result == 0);
        REQUIRE(compareFiles("../../test/chr20.malva.vcf", "out.vcf") == true);
    }
}