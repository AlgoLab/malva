#include <iostream>
#include "kseq.h"

#include "argument_parser.hpp"

#include "common.hpp"

int main(int argc, char *argv[]) {
    if (argc < 2)
    {
        cerr << "malva missing arguments" << endl;
        cerr << USAGE_MESSAGE << endl;
        return 1;
    }

    if (strncmp(argv[1], "index", 5) == 0)
    {
        return index_main(argc-1, argv+1);
    }
    else if (strncmp(argv[1], "call", 4) == 0)
    {
        return call_main(argc-1, argv+1);
    }
    else
    {
        cerr << "Could not interpret command " << argv[1] <<"." << endl;
        cerr << "Accepted commands are index and call." << endl;
        return 1;
    }
}