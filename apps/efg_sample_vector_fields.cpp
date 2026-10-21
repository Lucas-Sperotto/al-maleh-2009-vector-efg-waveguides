// SPDX-License-Identifier: MIT
// M6f: sample the ORIGINAL Cartesian EFG TM-H vector eigenfunctions at
// centroid coordinates exported by the independently compiled Whitney FEM.
// No modification, filtering, normalization or projection of K/M eigenpairs.
#include "al_maleh/cartesian.hpp"
#include "al_maleh/spectral.hpp"
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace al_maleh;
namespace {
struct Sample {
  std::size_t cell_id=0;
  double x=0,y=0;
};
std::vector<Sample> load_points(const std::filesystem::path& src){
  std::ifstream in(src);
  if(!in)throw std::runtime_error("Cannot read external FEM centroid CSV");
  std::string line;
  if(!std::getline(in,line)||line.rfind("cell_id,xc_m,yc_m,",0)!=0)
    throw std::runtime_error("Unexpected independent FEM centroid header");
  std::vector<Sample> samples;
  while(std::getline(in,line)){
    if(line.empty())continue;
    std::istringstream row(line);
    std::string a,b,c;
    if(!std::getline(row,a,',')||!std::getline(row,b,',')||
       !std::getline(row,c,','))
      throw std::runtime_error("Malformed FEM centroid CSV row");
    const auto idx=static_cast<std::size_t>(std::stoull(a));
    if(idx!=samples.size())
      throw std::runtime_error("FEM centroid row sequence not contiguous");
    // Native FEM edge_rect is 1x0.5; target EFG is 2x1. Both coordinates
    // must scale by 2, never just the spectral eigenvalue.
    const double x=2.0*std::stod(b),y=2.0*std::stod(c);
    if(!std::isfinite(x)||!std::isfinite(y)||
       x<0||x>2||y<0||y>1)
      throw std::runtime_error("FEM centroid outside rescaled guide");
    samples.push_back({idx,x,y});
  }
  if(samples.empty())throw std::runtime_error("No FEM centroid rows");
  return samples;
}
}
int main(int argc,char** argv){
  try{
    // nx ny native_TM_field.csv output.csv mode_index[6]
    if(argc!=11)throw std::invalid_argument(
      "Usage: efg_sample_vector_fields nx ny native_TM_field.csv output.csv idx1 idx2 idx3 idx4 idx5 idx6");
    const int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]);
    if(nx<4||ny<4||nx>17||ny>17||2*nx*ny>600)
      throw std::invalid_argument("Unsupported M6f EFG grid");
    std::array<int,6> indices{};
    for(int k=0;k<6;++k)indices[static_cast<std::size_t>(k)]=std::stoi(argv[5+k]);
    const auto samples=load_points(argv[3]);
    const Rectangle rectangle{0,2,0,1};
    const double radius=3.3*std::max(2.0/(nx-1),1.0/(ny-1));
    const MLSEvaluator mls(make_regular_cloud(rectangle,nx,ny,radius),1);
    const BackgroundQuadrature q{rectangle,nx-1,ny-1,3};
    const auto matrices=assemble_cartesian(mls,q);
    const auto spectrum=solve_generalized_symmetric(matrices.stiffness,matrices.mass);
    const auto n=static_cast<Eigen::Index>(mls.nodes().size());
    for(int k:indices)
      if(k<0||k>=spectrum.eigenvectors.cols())
        throw std::invalid_argument("M6f eigenvector index out of range");
    const std::filesystem::path output=argv[4];
    std::filesystem::create_directories(output.parent_path());
    std::ofstream out(output);
    if(!out)throw std::runtime_error("Cannot open M6f EFG sample output");
    out<<std::setprecision(17)<<"cell_id,x_efg,y_efg";
    for(int mode:indices)out<<",mode"<<mode<<"_hx,mode"<<mode<<"_hy";
    out<<"\n";
    for(const auto& sample:samples){
      const auto s=mls.evaluate({sample.x,sample.y});
      std::array<double,12> xy{};
      for(std::size_t k=0;k<s.active_nodes.size();++k){
        const auto node=static_cast<Eigen::Index>(s.active_nodes[k]);
        for(std::size_t m=0;m<6;++m){
          xy[2*m]+=s.phi[k]*spectrum.eigenvectors(node,indices[m]);
          xy[2*m+1]+=s.phi[k]*spectrum.eigenvectors(n+node,indices[m]);
        }
      }
      out<<sample.cell_id<<','<<sample.x<<','<<sample.y;
      for(double v:xy){
        if(!std::isfinite(v))throw std::runtime_error("Nonfinite reconstructed H");
        out<<','<<v;
      }
      out<<"\n";
    }
    if(!out)throw std::runtime_error("M6f EFG CSV write failed");
    std::cout<<"M6f EFG "<<nx<<"x"<<ny<<" EFG_nodes="<<n
             <<" FEM_centroids="<<samples.size()
             <<" selected_modes=6; raw spectral vectors unchanged\n";
  }catch(const std::exception& e){
    std::cerr<<"M6f ERROR: "<<e.what()<<"\n";return 1;
  }
}
