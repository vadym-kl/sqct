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


#ifndef SERIALIZERS_H
#define SERIALIZERS_H

#include <vector>
#include <fstream>

template< class T>
void read_stl( std::vector<T>& vec, const std::string& filename )
{
  std::ifstream ifs;
  ifs.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
  ifs.open( filename, std::ios_base::binary );

  std::size_t vec_sz = 0;
  ifs.read( (char*) &vec_sz, sizeof(std::size_t) );
  vec.resize(vec_sz);
  ifs.read( (char*) &vec.front(), sizeof(T) * vec_sz );
}

template< class T>
void write_stl( const std::vector<T>& s, const std::string& filename )
{
  if( s.size() == 0 )
    return;

  std::ofstream ofs;
  ofs.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
  ofs.open( filename, std::ios_base::binary );

  std::size_t set_sz = s.size();
  ofs.write( (const char*) &set_sz, sizeof(std::size_t) );
  ofs.write( (const char*) &( s.front() ), sizeof(T) * s.size() );
}

static std::vector<std::string> read_lines( const std::string& filename )
{
  std::ifstream ifs;
  ifs.exceptions ( std::ifstream::badbit );
  ifs.open(filename);

  std::vector<std::string> lines;
  while( !ifs.eof() )
  {
    std::string line;
    std::getline(ifs,line);
    if( line.size() != 0 && line.at(0) != '#' )
      lines.push_back(line);
  }
  return lines;
}

#endif // SERIALIZERS_H
