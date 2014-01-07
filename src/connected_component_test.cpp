#include <boost/config.hpp>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/timer/timer.hpp>
#include <boost/format.hpp>

int main(int argc, char *argv[])
{

  boost::timer::cpu_timer timer;

  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G;
  boost::add_edge(0, 6, G);
  boost::add_edge(1, 5, G);
  boost::add_edge(1, 8, G);
  boost::add_edge(2, 7, G);
  boost::add_edge(2, 9, G);
  boost::add_edge(3, 8, G);
  boost::add_edge(4, 6, G);
  boost::add_edge(4, 7, G);
  boost::add_edge(4, 9, G);

  std::vector<int> component(boost::num_vertices(G));
  int nConnectedComponents = boost::connected_components(G, &component[0]);

  boost::timer::cpu_times t = timer.elapsed();
  std::cout << boost::format("Elapsed time: %1% [ns]\n") % t.wall;

  printf("Number of components: %d\n", nConnectedComponents);
  for(int i = 0; i < component.size(); i++){
    printf("Vertex %d is in component %d\n", i, component[i]);
  }

  return 0;
}
