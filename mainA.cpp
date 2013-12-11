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
#include "appr/topt-bfs.h"
#include "test.h"

#include "requestprocessor.h"

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


////////////////////////////////////////////////////////////////////

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
    string file_name = "in.txt";
    try {

        po::options_description desc("Allowed options");
        desc.add_options()

            ("help,H", po::value< string >(&help_topic)->implicit_value(""),
             "Produce help message, see help <option name> for more details "
             "about specific option.")

//            ("test,T", po::value< string >(&(file_name))->implicit_value("generic"),
//             "Run tests")

            ("bfs,B", po::value< string >(&(file_name))->implicit_value("generic"),
             "Test that BFS and implemented algorithm produce the same results")

            ("gen,G", po::value< string >(&(file_name))->implicit_value("in.txt"),
             "Run unitaries approximation. Example: -G in ")

            ("about", "Information about the program.")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if( vm.count("about") ) {
            print_about_message();
            return 0;
        }

        if( vm.count("bfs") ) {
            bfs_results br;
            br.get();
            cup_test();
            return 0;
        }

//        if( vm.count("test") ) {
//            run_tests();
//            return 0;
//        }

        if( vm.count("help") ) {
            print_help(help_topic);
            return 0;
        }

        //by default we will start appropximating unitaries from in.txt
        {
          ifstream ifs(file_name);
          if( ifs )
          {
            std::vector<string> lines;
            while(ifs)
            {
              string line;
              getline(ifs,line);
              if(line.size() == 0)
                continue;
              if(line[0] == '#')
                continue;
              lines.push_back(line);
            }
            process_request(lines);
          }
          else
          {
            if( !print_help(help_topic) )
                cout << desc << endl;
            return 0;
          }
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

