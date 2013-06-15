//     Copyright (c) 2012 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com
//
//     This file is part of SQCT.
//
//     SQCT is free software: you can redistribute it and/or modify
//     it under the terms of the GNU Lesser General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     SQCT is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU Lesser General Public License for more details.
//
//     You should have received a copy of the GNU Lesser General Public License
//     along with SQCT.  If not, see <http://www.gnu.org/licenses/>.
//

#include "appr/toptzrot2.h"
#include "output.h"

#include "appr/zrot_cache.h"

#include <fstream>
#include <sstream>
#include <chrono>
#include <string>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <cmath>

// command line parcing code based on http://www.boost.org/doc/libs/1_49_0/doc/html/program_options/tutorial.html#id2499896
// and BOOST_ROOT/libs/program_options/example/first.cpp
namespace po = boost::program_options;
namespace btm = boost::timer;

using namespace std;

static bool print_help( const string& topic )
{
    map<string,string> hi;
    hi["in"];
    if( hi.count(topic) )
    {
        cout << hi[topic] << endl;
        return true;
    }

    return false;
}

////////////////////////////////////////////////////////////////////

static void print_about_message()
{
cout << "Copyright (c) 2013 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com" << endl << endl;
cout << "SQCT is free software: you can redistribute it and/or modify" << endl;
cout << "it under the terms of the GNU Lesser General Public License as published by" << endl;
cout << "the Free Software Foundation, either version 3 of the License, or" << endl;
cout << "(at your option) any later version." << endl;
cout << "" << endl;
cout << "SQCT is distributed in the hope that it will be useful," << endl;
cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
cout << "GNU Lesser General Public License for more details." << endl;
cout << "" << endl;
cout << "You should have received a copy of the GNU Lesser General Public License" << endl;
cout << "along with SQCT.  If not, see <http://www.gnu.org/licenses/>." << endl;
cout << "" << endl;
}

////////////////////////////////////////////////////////////////////

int main(int ac, char* av[])
{

    string help_topic;
    DiffApplicationParams diff_params;
    topt_app_params topt{80,"in","rz"};
    sqct_light_params sp;

    try {

        po::options_description desc("Allowed options");
        desc.add_options()

            ("help,H", po::value< string >(&help_topic)->implicit_value(""),
             "Produce help message, see help <option name> for more details "
             "about specific option.")

            ("gen,G", po::value< string >(&(topt.in_filename)),
             "File name with unitaries for approximation. Example: -G in ")

//            ("in,I", po::value< string >(&(topt.in_filename)),
//             "File name with unitaries for approximation. Example: -I sample")

            ("cache,C", po::value< string >(&(topt.out_filename))->default_value("rz.csv"),
             "File name with summary about "
             "all unitary approximation results.")

            ("out,O", po::value< string >(&(sp.out_filename))->default_value("out.qgl.xml"),
             "File name with "
             "unitary approximation results.")

            ("max-sde,M", po::value< int >(&(topt.max_sde))->default_value(80),
             "Maximal value of sde to use during approximation step.")

            ("diff,D", po::value< vector<string> >(&(diff_params.filenames) )->multitoken(),
             "Compares to caches of rotations approximations. "
             "Example: -D original/rz.csv rz.csv")

            ("about", "Information about the program.")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        sp.cache_filename = topt.out_filename;
        sp.in_filename = topt.in_filename;

        if( vm.count("diff") )
        {
          DiffApplication da(diff_params);
          da.run();
          return 0;
        }

        if( vm.count("about") ) {
            print_about_message();
            return 0;
        }

        if( vm.count("gen") ) {
            topt_app tapp(topt);
            tapp.run();
            return 0;
        }

        if( vm.count("in") ) {
            sqct_light_app lapp(sp);
            lapp.run();
            return 0;
        }

        {
            if( !print_help(help_topic) )
                cout << desc << endl;
            return 0;
        }

    }
    catch(exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!" << endl;
        return 1;
    }
}

