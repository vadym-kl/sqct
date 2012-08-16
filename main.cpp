//     Copyright (c) 2012 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com, Dmitri Maslov, Michele Mosca
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

#include "sk.h"
#include "exactdecomposer.h"
#include "gcommdecomposer.h"

#include "theoremverification.h"
#include "toptimalitytest.h"
#include "hoptimalitytest.h"

#include "netgenerator.h"

#include "output.h"

#include <fstream>
#include <sstream>
#include <chrono>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <cmath>

// command line parcing code based on http://www.boost.org/doc/libs/1_49_0/doc/html/program_options/tutorial.html#id2499896
// and BOOST_ROOT/libs/program_options/example/first.cpp
namespace po = boost::program_options;
namespace btm = boost::timer;

using namespace std;

//////////////////////////////////////////////////////////////////////////////

/// \brief Option of SK execution
struct SKOptions
{
    string in_filename; ///< input filename
    string out_filename; ///< output statistics filename
    bool math_fr; ///< should output be mathematica friendly
    int iterations_default; ///< number of SK iterations to perform by default
    bool write_dot_qc;///< should we output dot qc or not
    string out_dir;///< path where we should output dotqc files
    bool denom;///< how we interpret first element \b val of each input line ( Pi/\b val or just \b val)
};

//////////////////////////////////////////////////////////////////////////////

/// \brief Result of application of the SK algorithm
struct ApplicationResult
{
    circuit c;  ///< Cicrcuit
    int hc;     ///< Hadamard counts
    int tc;     ///< T and inverse of T gates counts
    int pc;     ///< Phase and inverse of Phase gate counts
    int plc;    ///< Pauli matrices counts
    int total;  ///< Total cost using cost function defined by vector gateLibrary::cost
    double dst; ///< Trace distance to approximation
    double tappr; ///< Time required to do an approximation
    double tdecomp; ///< Time required to do a decomposition
    int nr; ///< Number of iterations performed
    int denom_reduction; ///< Difference between \f$ sde(|\cdot|^2) \f$
                         /// before and after conversion to canonical form
    int denom ; ///< \f$ sde(|\cdot|^2) \f$ of resulting unitary

    /// \brief Summarize gate statistics for Clifford + T library
    void updateForCliffordT()
    {
        typedef gateLibrary gl;
        total = c.cost();
        auto counts = gl::toCliffordT( c.count() );
        hc = counts[ gl::H ];
        tc = counts[ gl::T ];
        pc = counts[ gl::P ];
        plc =  counts[ gl::Z ] + counts[ gl::X ] + counts[ gl::Y ];
    }

    /// \brief Outputs results in CSV format
    void toCSVStream( std::ostream& out ) const
    {
        out.precision(5);
        out.setf( ios_base::scientific );
        out << nr << ","
            << tc << ","
            << hc << ","
            << pc << ","
            << plc << ","
            << dst << ","
            << tappr << ","
            << tdecomp << "," << denom_reduction << "," << denom ;
    }

    /// \brief Outputs results in Mathematica friendly format
    void toMathematicaStream( std::ostream& out ) const
    {
        out.precision(10);
        out.setf( ios_base::fixed ); // use fixed as mathematica does not understand scientific c++ format
        out << "AppResult["
            << nr << ","
            << tc << ","
            << hc << ","
            << pc << ","
            << plc << ",";
        out.precision(60);
        out << dst << "," ;
        out.precision(10);
        out << tappr << ","
            << tdecomp << "," << denom_reduction << "," << denom << "]";
    }

    /// \brief Outputs dot qc file with circuit and adds information about genration process into comments
    void generateDotQc( std::string filename, std::string circuit_name, std::string symbolic_form )
    {
        ofstream ofs(filename.c_str());
        ofs << "# Copyright (c) 2012 Vadym Kliuchnikov, Dmitri Maslov, Michele Mosca" << endl;
        ofs << "# Computed by QSQ, based on arXiv:1206.5236" << endl;
        ofs << "# Published at: http://qcirc.iqc.uwaterloo.ca" << endl;
        //ofs << "# Software published at: http://qcirc.iqc.uwaterloo.ca" << endl;
        ofs << "# Symbolic form of unitary to approximate:" << symbolic_form << endl;
        ofs << "# Trace distance between unitary and approximation:" << dst << endl;
        ofs << "# Total number of gates:" << tc + hc + pc + plc << endl;
        ofs << "# T and T^{Dagger} gates:" << tc << endl;
        ofs << "# Hadamard gates:" << hc << endl;
        ofs << "# Phase and Phase^{Dagger} gates:" << pc << endl;
        ofs << "# Pauli gates:" << plc << endl;
        ofs << "# Total cost: " << total << endl;
        ofs << "# Number of iterations:" << nr << endl;
        ofs << "# Approximation time(seconds):" << tappr << endl;
        ofs << "# Decomposition time(seconds):" << tdecomp << endl;
        ofs << "# Reduction:" << denom_reduction << endl;
        ofs << "# sde(abs(z)^2):" << denom << endl;

        ofs << ".v 1" << endl;
        ofs << ".i 1" << endl;
        ofs << ".o 1" << endl;
        ofs << "BEGIN " << circuit_name <<"(1)" << endl;
        c.toStream(ofs);
        ofs << "END " << circuit_name << endl;
        ofs << "BEGIN" << endl;
        ofs << circuit_name << " 1" << endl;
        ofs << "END" << endl;
        ofs.close();
    }

};

///////////////////////////////////////////////////////////////////

/// \brief Process input file with inputs for SK in batch mode, outputs statistics etc
class SKApplication
{
public:
    SKApplication( SKOptions& skopt ) :
        opt( skopt )
    {}

    /// \brief Finds circuit and all supplementary information for unitary
    void calculate( const matrix2x2hpr& matrix, int recursion_level, ApplicationResult& ar )
    {
        boost::timer::cpu_times ct;
        boost::timer::cpu_times ct2;

        sk::Me res;
        {
            boost::timer::auto_cpu_timer t;
            skd.decompose(matrix,res,recursion_level);
            ct = t.elapsed();
        }

        int gde_before = res.min_gde_abs2();
        res.reduce();
        int gde_after = res.min_gde_abs2();
        ar.denom_reduction = gde_before - gde_after;
        ar.denom = res.max_sde_abs2();

        {
            boost::timer::auto_cpu_timer t;
            ed.decompose(res,ar.c);
            ct2 = t.elapsed();
        }

        sk::Ma conv(res);

        ar.dst = trace_dist(matrix,conv);
        ar.tappr = (double) ct.wall * 1e-9;
        ar.tdecomp = (double) ct2.wall * 1e-9;
        ar.nr = recursion_level;
        ar.updateForCliffordT();
    }

    /// \brief Process input file defined by opt taking into account specified options
    void process()
    {
        static const double TwoPi = 2. * hprHelpers::toMachine( hprHelpers::pi() );
        ifstream ifs( opt.in_filename.c_str() );
        ofstream ofs( opt.out_filename.c_str() );
        Rotation r;
        int iter = 0;
        ApplicationResult ar;
        int line_number = 0;
        while( true )
        {
            string line;
            getline(ifs,line);
            iter = opt.iterations_default;
            istringstream ss(line,istringstream::in);
            if( ! ifs ) break;

            double val;
            ss >> val >> r.nx >> r.ny >> r.nz >> iter;
            if( opt.denom )
            {
                r.num = 1.0;
                r.den = val;
            }
            else
            {
                r.num = val;
                r.den = TwoPi;
            }

            calculate( r.matrix() , iter, ar );

            stringstream fname;
            fname <<  opt.out_dir << "/";
            if( r.isSpecial() != 'N' )
                fname << r.name() << "-" << iter;
            else
                fname << opt.out_filename << "." << line_number;
            fname << ".qc" ;

            if( opt.write_dot_qc )
                ar.generateDotQc( fname.str() ,
                               r.name(),
                               r.symbolic() );

            if( opt.math_fr )
            {
                ofs << "{" << r.Mathematica() << ",";
                ar.toMathematicaStream(ofs);
                ofs << "}" << endl;
            }
            else
            {
                ofs << r.CSV() << ",";
                ar.toCSVStream(ofs);
                ofs << endl;
            }

            line_number++;
        }
    }
private:
    const SKOptions& opt;
    sk skd;
    exactDecomposer ed;
};

///////////////////////////////////////////////////////////////////

bool print_help( const string& topic )
{
    map<string,string> hi;
    if( hi.count(topic) )
    {
        cout << hi[topic] << endl;
        return true;
    }

    return false;
}

/////////////////////////////////////////////////////////////////////

struct ResynthOptions
{
    vector<string> circuit_files;
    bool reverse;
    bool math_fr;
};

////////////////////////////////////////////////////////////////////

/// \brief Resyntheisizes circuits using exact decomposition algoritm
struct ResynthApplication
{
    ResynthApplication( const ResynthOptions& opt ) :
        ro(opt)
    {}

    void process()
    {
        for( const auto& fn : ro.circuit_files )
            process(fn);
    }

    void process( const string& filename )
    {
        const gateLibrary& gl = gateLibrary::instance();
        ifstream ifs(filename.c_str());
        ofstream ofs( (filename + ".out").c_str() );
        // read circuit first
        circuit in, out;
        in.fromStream(ifs,ro.reverse);
        ed.decompose(in,out);

        auto in_cost = in.count();
        int in_total = in.cost();
        in_cost = gl.toCliffordT(in_cost);

        auto out_cost = out.count();
        int out_total = out.cost();
        out_cost = gl.toCliffordT(out_cost);

        string cb = "# ";
        string ce = "";
        if( ro.math_fr )
            cb = "(* ",ce = " *)";

        ofs << cb << "Computed by ESK, based on arXiv:1206.5236" << ce << endl;
        ofs << cb << "Software published at: http://qcirc.iqc.uwaterloo.ca"<< ce << endl;
        ofs << cb << "Gate counts for input circuit [T+Td,H,P+Pd,X,Z,Y] : ["
            << in_cost[gl.T] << "," << in_cost[gl.H] << ","
            << in_cost[gl.P] << "," << in_cost[gl.X] << ","
            << in_cost[gl.Z] << "," << in_cost[gl.Y] << "]" << ce << endl;
        ofs << cb << "Total cost of input:" << in_total << ce << endl;

        ofs << cb << "Gate counts for output circuit [T+Td,H,P+Pd,X,Z,Y] : ["
            << out_cost[gl.T] << "," << out_cost[gl.H] << ","
            << out_cost[gl.P] << "," << out_cost[gl.X] << ","
            << out_cost[gl.Z] << "," << out_cost[gl.Y] << "]" << ce << endl;
        ofs << cb << "Total cost of output:" << out_total << ce << endl;

        if( ro.math_fr )
            out.toMathStream(ofs);
        else
            out.toStreamSym(ofs);

        ofs.close();
        ifs.close();
    }

    exactDecomposer ed;
    const ResynthOptions& ro;
};

////////////////////////////////////////////////////////////////////

/// \brief Options for epsilon net generation
struct enetOptions
{
    /// \brief If one element specified -- upper bound for sde of epsilon net to be
    /// generated. If two elements specified -- interval of sde to be generated.
    vector<int> epsilon_net_layers;
};

////////////////////////////////////////////////////////////////////

struct enetApplication
{
    enetApplication( const enetOptions& options ) :
        m_options(options)
    {}

    /// \brief Checks if all layers generated on initial state are available
    void check_initial()
    {
        initial_ok = true;
        if( m_layers[0] == 0 )
            initial_ok = false;

        for( int i = 2; (i < initial_end) && initial_ok ; ++i )
            if( m_layers[i] == 0 ) initial_ok = false;
    }

    /// \brief Collects information about parts of epsilon net available
    void check_files()
    {
         m_layers.resize(100,0);
         // collect information about available layers
         for( int i = 0; i < 100; ++i )
         {
             string name = netGenerator::fileName(i);
             ifstream ifs(name);
             if( !ifs )
                 m_layers[i] = 0;
             else
                 m_layers[i] = 1;
         }
    }

    /// \brief Generates all layers with \f$ sde(|\cdot|^2) \f$ in the interval [st,end)
    void generate( int st, int end )
    {
        if( ! initial_ok && st < initial_end )
            m_ng.generateInitial();

        for( int i = max(initial_end,st) ; i < end; ++i )
        {
            if( m_layers[i] == 0 )
            {
                epsilonnet base_net;
                base_net.loadFromFile( netGenerator::fileName(i-1).c_str() );
                unique_ptr<epsilonnet> res( netGenerator::generate( base_net ) );
                res->saveToFile( netGenerator::fileName(i).c_str() );
            }
        }
    }

    /// \brief Do the work
    void process()
    {
        check_files();
        check_initial();

        int sz  = m_options.epsilon_net_layers.size();
        int b = m_options.epsilon_net_layers[0];
        switch( sz )
        {
        case 1:
            generate(0,b);
            break;
        case 2:
            generate(b,m_options.epsilon_net_layers[1]);
            break;
        default:
            std::cerr << "Wrong number of parameters. See help." << endl;
        }
    }

    const enetOptions& m_options;       ///< Application options
    std::vector<int> m_layers;          ///< Ones for available layers, zeros for not availible layers
    static const int initial_end;  ///< Maximal layer generated on initial state
    bool initial_ok;                    ///< If all initial layers are availible
    netGenerator m_ng;                  ///< Epsilon net generator
};

const int enetApplication::initial_end = 21;

////////////////////////////////////////////////////////////////////

struct RandomRotationsOpts
{
    string filename;
    int samples;
};

////////////////////////////////////////////////////////////////////

struct RotationGeneratorApp
{
    RotationGeneratorApp( RandomRotationsOpts& opts) :
        m_opt(opts)
    {}
    /// \brief Generates file with rotations by random angle around random axis
    void process()
    {
        ofstream ofs( m_opt.filename.c_str() );
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        static const double PI = 3.1415926535897932384626433832795028841971693993751;
        normal_distribution<double> nd(0.0,1.0);
        uniform_real_distribution<double> urd( 0.0, 2*PI );

        for( int i = 0; i < m_opt.samples; ++i )
        {
            double x = nd(generator);
            double y = nd(generator);
            double z = nd(generator);
            double norm = sqrt(x*x + y*y + z*z);
            x /= norm;
            y /= norm;
            z /= norm;
            double phi = urd(generator);
            ofs << phi << " " << x << " " << y << " " << z << endl;
        }
    }

    const RandomRotationsOpts& m_opt;
};

////////////////////////////////////////////////////////////////////

void print_about_message()
{
cout << "Copyright (c) 2012 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com" << endl << endl; 
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
cout << "The program code based on results of http://arxiv.org/abs/1206.5236." << endl;
cout << "It also implements the version of Solovay Kitaev algorithm described" << endl;
cout << "in http://arxiv.org/abs/quant-ph/0505030. " << endl;
cout << "" << endl; 
cout << "The list of used libraries, corresponding licenses and authors:" << endl<< endl;
cout << "Boost | Boost Software License" << endl<< endl;
cout << "The GNU Multiple Precision Arithmetic Library | GNU Lesser General Public License v3" << endl;
cout << "(by Free Software Foundation, Inc.)" << endl<< endl;
cout << "The GNU MPFR Library the library | GNU Lesser General Public License v3" << endl;
cout << "(by Free Software Foundation, Inc.," << endl;
cout << " Contributed by the AriC and Caramel projects, INRIA. )" << endl << endl;
cout << "mpfr::real | GNU Lesser General Public License v3" << endl;
cout << "(by Christian Schneider <software(at)chschneider(dot)eu> )" << endl;
cout << "" << endl;
cout << "Source code of this program can be obtained from: " << endl;
}

////////////////////////////////////////////////////////////////////

int main(int ac, char* av[])
{
    // help
    string help_topic = "";
    // sk application parameters
    SKOptions sko;
    // epsilon net generation
    enetOptions eopt;
    // resynthesise
    ResynthOptions ro;
    // randomized rotations
    RandomRotationsOpts rro;

    try {

        po::options_description desc("Allowed options");
        desc.add_options()

            ("help,H", po::value< string >(&help_topic)->implicit_value(""),
             "Produce help message, see help <option name> for more details "
             "about specific option.")

            ("in,I", po::value< string >(&(sko.in_filename)),
             "File name with unitaries for approximation.")

            ("iterations,N", po::value< int >(&(sko.iterations_default))->default_value(3),
             "Default number of iteration used by the Solovay Kitaev algorithm, "
             "when not specified in input file.")

            ("denominator,R", po::value< bool >(&(sko.denom))->default_value(true),
             "When true program interprets first entry in each line of input file as "
             "denominator of Pi, and as rotation angle otherwise.")

            ("out,O", po::value< string >(&(sko.out_filename))->default_value("out.csv"),
             "File name with summary about "
             "unitary approximation results.")

            ("math-fr,F", po::value< bool >(&(sko.math_fr))->default_value(false),
             "Should output be Mathematica friendly. "
             "Default output format is CSV.")

            ("dotqc,C", po::value< bool >(&(sko.write_dot_qc))->default_value(false),
             "Should we output dotqc file for each generated circuit.")

            ("out-dir,D", po::value< string >(&(sko.out_dir))->default_value("out"),
             "Directory to output circuits.")

            ("epsilon-net,E", po::value< vector<int> >(&(eopt.epsilon_net_layers))->multitoken(),
             "Generates epsilon net layers. If one values specified generates all layers "
             "with sde(|.|^2) less than this value. If two values specified then sde(|.|^2) "
             "of result lies in interval [fist,second) ")

            ("theory-topt",
             "Verifies conjecture about T optimality ")

            ("theory-hopt",
             "Performs brute force check necessary for proof of result in Appendix 2 "
             "in http://arxiv.org/abs/1206.5236")

            ("theory-correctness",
             "Performs brute force check described by Algoritm 2 "
             "in http://arxiv.org/abs/1206.5236 ( proves Lemma 3 ).")

            ("resynth,S", po::value< vector<string> >(&(ro.circuit_files) )->multitoken(),
             "Resynthesize circuit. Output may differ from original by global phase.")

            ("reverse", po::value< bool >(& (ro.reverse) )->default_value(false),
             "True if circuit file should be interpreted in matrix multiplication order.")

            ("random-rot", po::value< int >(&(rro.samples)),
             "Generates file with specified number of random rotations.")

            ("about", "Information about the program.")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if( vm.count("resynth") ) {
            ro.math_fr = sko.math_fr;
            ResynthApplication ra(ro);
            ra.process();
            return 0;
        }

        if(vm.count("epsilon-net"))
        {
            enetApplication eapp(eopt);
            eapp.process();
            return 0;
        }

        if( vm.count("theory-correctness") ) {
            cout << is_theorem_true() << endl;
            return 0;
        }

        if( vm.count("theory-topt") ) {
            toptimalitytest tt;
            return 0;
        }

        if( vm.count("theory-hopt") ) {
            hoptimalitytest ht;
            return 0;
        }

        if ( vm.count("random-rot") ) {
            rro.filename = sko.out_filename;
            RotationGeneratorApp rga(rro);
            rga.process();
            return 0;
        }

        if( vm.count("about") ) {
            print_about_message();
            return 0;
        }

        if ( vm.count("in") ) {
            SKApplication app(sko);
            app.process();
            return 0;
        }

        if ( true ) {
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
