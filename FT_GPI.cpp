/* FT-GPI
 * Copyright (C) 2022 Rodrigo Canovas
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/ .
 * */


#include <unistd.h>
#include <iostream>
#include "./include/ftgpi.h"
#include <boost/program_options.hpp>

int main(int argc, char* argv[]) {

    std::string file, out_file_pred;
    std::string training_file = "";
    double d = 1000, e = 1000, k = 1.0;
    int n = 1000, c = 1800, i = -1;
    int m = 17, M = 35;
    int t1 = 80, t2 = 200, t3 = 500, t4 = 80;
    bool as_table = false, d_type = false, both_models = false;
    namespace po = boost::program_options;
    po::options_description desc("FT-GPI Options");
    desc.add_options()
            // First parameter describes option name/short name
            // The second is parameter to option
            // The third is description
            ("help,h", "print usage message")
            ("output,o", po::value<std::string>(), "pathname for output. Default: file_name.gpi")
            ("training,t", po::value<std::string>(), "pathname of the training set of GPI proteins. Default: none")
            ("as-table,u", "if set, prints results as a table rather than as a list. Default: not set")
            ("as-percentages,p", "if set, the start and end parameters are considered as percentages over the size of the sequence. Default: not set -> fixed integer distances")
            ("both,b", "if set, checks proteins alternative model if strong model fails. Default: program only checks strong models")
            ("start,d", po::value<double>(), "maximum distance between the Nterm with the start of the sequence. Default: 15 amino acids (or 5% if p is set)")
            ("end,e", po::value<double>(), "maximum distance between the Cterm with the end of the sequence. Default: 10 amino acids (or 1% if p is set)")
            ("n-score,n", po::value<int>(), "minimum score of the Nterm. Default: 1000")
            ("c-score,c", po::value<int>(), "minimum score of the Cterm. Default: 1800")
						("i-score,i", po::value<int>(), "maximum score allowed for all internal helixes. Default: -1 (not seted)")
            ("i-score-percentage,k", po::value<double>(), "maximum (excluded) allowed score for internal helixes computed as k*max-score(Nterm, Cterm). Default: 1.0");
    po::options_description tmpred_op("TMPred Options");
    tmpred_op.add_options()
            ("min-len,m", po::value<int>(), "minimal length of transmembrane sequence. Default: 17")
            ("max-len,M", po::value<int>(), "maximal length of transmembrane sequence. Default: 35")
            ("low-osl,l", po::value<int>(), "low orientational significance level. Default: 80")
            ("high-osl,g", po::value<int>(), "high orientational significance level. Default: 200")
            ("tm-osl,s", po::value<int>(), "TM-existence significance level. Default: 500")
            ("avg-osl,a", po::value<int>(), "average orientation significance level. Default: 80");
    desc.add(tmpred_op);
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help") or argc < 2) {
        std::cout << "Usage: " << argv[0] << " file_name <opt>" << std::endl;
        std::cout << desc << std::endl;
        return 0;
    }
    file = argv[1];
    out_file_pred = file;

    if (vm.count("output")) { out_file_pred = vm["output"].as<std::string>();}
    if (vm.count("training")) { training_file = vm["training"].as<std::string>();}
    if (vm.count("as-table")) { as_table = true; }
    if (vm.count("as-percentages")) { d_type = true; }
    if (vm.count("both")) { both_models = true; }
    if (vm.count("start")) { d = vm["start"].as<double>(); }
    if (vm.count("end")) { e = vm["end"].as<double>(); }
    if (vm.count("n-score")) { n = vm["n-score"].as<int>(); }
    if (vm.count("c-score")) { c = vm["c-score"].as<int>(); }
		if (vm.count("i-score")) { i = vm["i-score"].as<int>(); }
    if (vm.count("i-score-percentage")) { k = vm["i-score-percentage"].as<double>(); }
    if (vm.count("min-len")) { m = vm["min-len"].as<int>();}
    if (vm.count("max-len")) { M = vm["max-len"].as<int>();}
    if (vm.count("low-osl")) { t1 = vm["low-osl"].as<int>();}
    if (vm.count("high-osl")) { t2 = vm["high-osl"].as<int>();}
    if (vm.count("tm-osl")) { t3 = vm["tm-osl"].as<int>();}
    if (vm.count("avg-osl")) { t4 = vm["avg-osl"].as<int>();}


    if (d_type) {
        if (d == 1000) d = 5;
        if (e == 1000) e = 1;
        fdgpi::ftgpi<double>(file, out_file_pred, training_file,
                               m, M, t1, t2, t3, t4, d, e, n, c, i, k, d_type, both_models, as_table);
    }
    else {
        if (d == 1000) d = 15;
        if (e == 1000) e = 10;
        fdgpi::ftgpi<int>(file, out_file_pred, training_file,
                            m, M, t1, t2, t3, t4, d, e, n, c, i, k, d_type, both_models, as_table);
    }

    return 0;
}

