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

// VALUE of TOLERANCE RANGE MATCH of GQ
int tolerance = 0;

/**
 *1. Check if the arguments are at least three
 *2. If 3 execute compare_vcf
 *3. If 4 throw an error
 *4. If 5 check if the Option is "-t", else throw an error
 *5. If the Option value is accepted execute compare_vcf, else throw an error
 **/
int main(int argc, char *argv[])
{
    // We expect at least three arguments
    if (argc < 3)
    {
        std::cerr << "<< MISSING ARGUMENTS Mandatory: <output_malva.vcf> <sample.vcf> >>" << std::endl;
        return 1;
    }
    else if (argc == 3)
    {
        //[0]Program Name, [1]Geno_path, [2]Sample_path
        compare_vcf(argv[1], argv[2]);
    }
    else if (argc == 4)
    {
        std::cerr << "<< MISSING ARGUMENTS: -t [t_value] <output_malva.vcf> <sample.vcf> >>" << std::endl;
        return 1;
    }
    else if (argc == 5)
    {
        //[0]Program Name, [1]Option, [2]Opt_Value, [3]Geno_path, [4]Sample_path
        if (strcmp(argv[1], "-t") != 0)
        {
            std::cerr << "<< WRONG OPTION COMMAND, try -t >>" << std::endl;
            return 1;
        }
        // tolerance min: 0 (default value), max: 100
        tolerance = atoi(argv[2]);
        if ((tolerance < 0) || (tolerance > 100))
        {
            std::cerr << "<< WRONG OPTION VALUE: Must be between 0 and 100 >>" << std::endl;
            return 1;
        }
        compare_vcf(argv[3], argv[4]);
    }
    else
    {
        return 1;
    }
    return 0;
}