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
