#include <boost/config.hpp>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/timer/timer.hpp>
#include <boost/format.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

int main(int argc, char *argv[])
{

  boost::mt19937 rng;
  boost::uniform_01<> ud;
  boost::variate_generator<boost::mt19937, boost::uniform_01<> > gen(rng, ud);

  int const C_size = 5;
  double** C = new double*[C_size];
  for(int i = 0; i < C_size; i++){
    C[i] = new double[C_size];
    for(int j = 0; j < C_size; j++){
      C[i][j] = 0;
    }
  }
  C[0][3] = gen();
  C[1][0] = gen();
  C[2][2] = gen();
  C[2][4] = gen();
  C[4][1] = gen();
  C[4][2] = gen();
  C[4][4] = gen();

  for(int i = 0; i < C_size; i++){
    for(int j = 0; j < C_size; j++){
      printf("%f  ", C[i][j]);
    }
    printf("\n");
  }


  boost::timer::cpu_timer timer;

  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G;
  for(int i = 0; i < C_size; i++){
    for(int j = 0; j < C_size; j++){
      if (C[i][j] != 0 ){
	boost::add_edge(i, j+C_size, G);
      }
    }
  }
  std::vector<int> cc_results(boost::num_vertices(G));
  int nConnectedComponents = boost::connected_components(G, &cc_results[0]);

  std::vector<unsigned int> components_m [nConnectedComponents];
  for(int i = 0; i <  C_size; i++){
    components_m[cc_results[i]].push_back(i);
  }
  std::vector<unsigned int> components_z [nConnectedComponents];
  for(int i = C_size; i < C_size*2; i++){
    components_z[cc_results[i]].push_back(i - C_size);
  }

  boost::timer::cpu_times t = timer.elapsed();
  std::cout << boost::format("Elapsed time: %1% [ns]\n") % t.wall;

  for(int i = 0; i < nConnectedComponents; i++){
    printf("Component %d: ", i);
    for(int j = 0; j < components_m[i].size(); j++){
	printf("m%d ", components_m[i][j]);
    }
    for(int j = 0; j < components_z[i].size(); j++){
	printf("z%d ", components_z[i][j]);
    }
    printf("\n");
    for(int m = 0; m < components_m[i].size(); m++){
      for(int n = 0; n < components_z[i].size(); n++){
	printf("%f  ", C[ components_m[i][m] ][ components_z[i][n] ]);
      }  
      printf("\n");
    }
  }

  for(int i = 0; i < C_size; i++){
    delete[] C[i];
  }
  delete[] C;

  return 0;
}
