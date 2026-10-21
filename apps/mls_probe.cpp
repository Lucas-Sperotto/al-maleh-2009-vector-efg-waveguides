#include "al_maleh/mls.hpp"
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
  std::vector<al_maleh::Node> cloud;
  for (int j = -2; j <= 2; ++j)
    for (int i = -2; i <= 2; ++i)
      cloud.push_back({{i*0.5, j*0.5}, 2.5});
  const al_maleh::MLSEvaluator e(cloud, 2);
  const auto r = e.evaluate({0.2, -0.15});
  std::cout << std::setprecision(16);
  std::cout << "# polynomial_degree=2; weight=wendland_c2; condition_A="
            << r.condition_number << '\n';
  std::cout << "node_index,phi,dx,dy,dxx,dxy,dyy\n";
  for (std::size_t i=0; i<r.active_nodes.size(); ++i)
    std::cout << r.active_nodes[i] << ',' << r.phi[i] << ','
              << r.dx[i] << ',' << r.dy[i] << ',' << r.dxx[i] << ','
              << r.dxy[i] << ',' << r.dyy[i] << '\n';
}
