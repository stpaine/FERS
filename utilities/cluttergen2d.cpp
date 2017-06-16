/// Two dimensional clutter generator, this is an extremely over-simplified program, but it does the job

#include <iostream>
#include <fstream>
#include <boost/random.hpp>

using namespace std;

int main()
{
  int samples;
  cout << "Number of clutter samples: ";
  cin >> samples;
  double start_range_x;
  cout << "Start range x: ";
  cin >> start_range_x;
  double range_x;
  cout << "Range x: ";
  cin >> range_x;

  double start_range_y;
  cout << "Start range y: ";
  cin >> start_range_y;
  double range_y;
  cout << "Range y: ";
  cin >> range_y;

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
  boost::uniform_real<double> ud_x(start_range_x, start_range_x+range_x);
  boost::uniform_real<double> ud_y(start_range_y, start_range_y+range_y);
  boost::normal_distribution<double> nd(0, spread);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > gen_x(rng, ud_x);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > gen_y(rng, ud_y);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > sprgen(rng, nd);
  fo << "<incblock>";
  for (int i = 0; i < samples; i++)
    {
      fo << "<platform name=\"clutter\">\n";
      fo << "<motionpath interpolation=\"cubic\">\n";
      double pos_x = gen_x();
      double pos_y = gen_y();
      fo << "<positionwaypoint>\n<x>" << pos_x << "</x>\n<y>" << pos_y << "</y>\n<altitude>0</altitude>\n<time>0</time>\n</positionwaypoint>\n";
      fo << "<positionwaypoint>\n<x>" << pos_x+time*sprgen() << "</x>\n<y>" << pos_y+time*sprgen() << "</y>\n<altitude>0</altitude>\n<time>" << time << "</time>\n</positionwaypoint>\n";
      fo << "</motionpath>\n";
      fo << "<fixedrotation><startazimuth>0.0</startazimuth><startelevation>0.0</startelevation><azimuthrate>0</azimuthrate><elevationrate>0</elevationrate></fixedrotation>\n";
      fo << "<target name=\"wings\">\n<rcs type=\"isotropic\">\n<value>";
      fo << rcs;
      fo << "</value>\n</rcs>\n</target>\n</platform>\n\n";
    }
  fo << "</incblock>";
}
