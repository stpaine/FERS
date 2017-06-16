/// One dimensional clutter generator, this is an extremely over-simplified program

#include <iostream>
#include <fstream>
#include <boost/random.hpp>

using namespace std;

int main()
{
  int samples;
  cout << "Number of clutter samples: ";
  cin >> samples;
  double start_range;
  cout << "Start range: ";
  cin >> start_range;
  double range;
  cout << "Range: ";
  cin >> range;
  double rcs;
  cout << "RCS: ";
  cin >> rcs;
  double spread;
  cout << "Stdev of spreading: ";
  cin >> spread;
  double time = 0;
  if (spread != 0) {
    cout << "Simulation end time";
    cin >> time;
  }
  string filename;
  cout << "Filename:";
  cin >> filename;
  //Open the file
  ofstream fo(filename.c_str());
  //Create the rng
  boost::mt19937 rng;
  boost::uniform_real<double> ud(start_range, start_range+range);
  boost::normal_distribution<double> nd(0, spread);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > gen(rng, ud);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > sprgen(rng, nd);
  fo << "<incblock>";
  for (int i = 0; i < samples; i++)
    {
      fo << "<platform name=\"clutter\">\n";
      fo << "<motionpath interpolation=\"cubic\">\n";
      double pos = gen();
      fo << "<positionwaypoint>\n<x>" << pos << "</x>\n<y>0</y>\n<altitude>0</altitude>\n<time>0</time>\n</positionwaypoint>\n";
      fo << "<positionwaypoint>\n<x>" << pos+time*sprgen() << "</x>\n<y>0</y>\n<altitude>0</altitude>\n<time>" << time << "</time>\n</positionwaypoint>\n";
      fo << "</motionpath>\n";
      fo << "<fixedrotation><startazimuth>0.0</startazimuth><startelevation>0.0</startelevation><azimuthrate>0</azimuthrate><elevationrate>0</elevationrate></fixedrotation>\n";
      fo << "<target name=\"wings\">\n<rcs type=\"isotropic\">\n<value>";
      fo << rcs;
      fo << "</value>\n</rcs>\n</target>\n</platform>\n\n";
    }
  fo << "</incblock>";
}
