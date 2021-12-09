/**
 * MALVA-TEST is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MALVA-TEST is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MALVA-TEST; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include "malva_test.hpp"

int main(int argc, char *argv[]) {
    //We expect 3 arguments: the program name, the geno path and the sample path
    if (argc < 3) {
        std::cerr << "MISSING ARGUMENTS" << std::endl;
        return 1;
    }
    compare_vcf(argv[1],argv[2]);
    return 0;
}