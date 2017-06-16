#include <iostream>
#include "rspath.h"
#include "cycle.h"

using namespace std;
using namespace rs;

int main() {
  ticks t1, t2;
  
  Path testpath(Path::RS_INTERP_CUBIC);
  Coord coord = {0, 0, 0, 0};
  //Add
  testpath.AddCoord(coord);
  coord.x = 1; coord.t = 1;
  testpath.AddCoord(coord);
  coord.x = 3; coord.t = 2;
  testpath.AddCoord(coord);
  coord.x = -1; coord.t = 3;
  testpath.AddCoord(coord);

  t1 = getticks();
  testpath.Finalize();
  t2 = getticks();
  cout << "Setup took " << elapsed(t2, t1) << " ticks" << endl;

  t1 = getticks();
  for (int i = -100; i < 3100; i++)
    testpath.GetPosition(i/1000.0, coord);
  t2 = getticks();
  cout << "Interp took " << elapsed(t2, t1) << " ticks" << endl;

  for (int i = 0; i <= 30; i+=2) {
    testpath.GetPosition(i/10.0, coord);
    cout << coord.t << " " << coord.x << " " << coord.y << " " << coord.z << endl;
  }
  return 0;
}
